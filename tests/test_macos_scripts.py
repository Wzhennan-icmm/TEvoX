#!/usr/bin/env python3
"""Test Mac script orchestration with fake Apple tools, never Mach-O execution.

These tests use real shell, file, archive and checksum operations. Compiler,
SDK, architecture and signing tools are simulated; native Mac CI is required
for Apple Clang, SDK, signature and executable validation.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
MOCK_TOOL = r'''
import json, os, pathlib, sys
name = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
with open(os.environ['MOCK_TOOL_LOG'], 'a') as log:
    log.write(json.dumps([name] + args) + '\n')

def emit_binary(path, architectures):
    path.write_text('#!/bin/sh\n# mock_archs: ' + ' '.join(architectures) +
                    '\nprintf "Usage: mock TEvoX\\n"\n')
    path.chmod(0o755)

def archs(path):
    text = pathlib.Path(path).read_text()
    return next(line.split(': ', 1)[1].split() for line in text.splitlines()
                if line.startswith('# mock_archs: '))

if name == 'uname':
    print(os.environ.get('MOCK_OS', 'Darwin') if args == ['-s'] else
          os.environ.get('MOCK_ARCH', 'arm64'))
elif name == 'sw_vers':
    print(os.environ.get('MOCK_OS_VERSION', '15.0'))
elif name == 'sysctl':
    print(os.environ.get('MOCK_ROSETTA', '0'))
elif name == 'xcrun':
    if '--show-sdk-path' in args:
        print(os.environ['MOCK_SDK'])
    elif '--find' in args:
        print(pathlib.Path(sys.argv[0]).parent / args[args.index('--find') + 1])
    else:
        sys.exit('Unsupported mock xcrun invocation')
elif name == 'clang':
    if args == ['--version']:
        print('Apple clang mock (NOT a real Apple compiler)')
    else:
        if os.environ.get('MOCK_COMPILE_FAIL') == '1':
            sys.exit(7)
        architecture = args[args.index('-arch') + 1]
        output = pathlib.Path(args[args.index('-o') + 1])
        assert '-mmacosx-version-min=11.0' in args
        sdk = args[args.index('-isysroot') + 1]
        assert pathlib.Path(sdk).is_dir()
        output.parent.mkdir(parents=True, exist_ok=True)
        if '-c' in args:
            source = pathlib.Path(args[args.index('-c') + 1])
            assert source.is_file()
            output.write_text(json.dumps({'arch': architecture}))
        else:
            for arg in args:
                if arg.endswith('.o'):
                    assert json.loads(pathlib.Path(arg).read_text())['arch'] == architecture
            emit_binary(output, [architecture])
elif name == 'lipo':
    if '-create' in args:
        merged = []
        for path in args[1:args.index('-output')]:
            merged.extend(archs(path))
        emit_binary(pathlib.Path(args[args.index('-output') + 1]), merged)
    elif '-verify_arch' in args:
        assert set(args[1:-1]).issubset(archs(args[-1]))
    elif '-archs' in args:
        print(' '.join(archs(args[-1])))
    else:
        sys.exit('Unsupported mock lipo invocation')
elif name == 'codesign':
    binary = pathlib.Path(args[-1])
    if os.environ.get('MOCK_SIGN_FAIL') == '1':
        sys.exit(8)
    if '--verify' in args:
        assert '# mock signed' in binary.read_text()
    else:
        assert args[args.index('--sign') + 1] == '-'
        with binary.open('a') as out:
            out.write('# mock signed\n')
else:
    sys.exit('Unknown simulated tool: ' + name)
'''


class MacScriptOrchestrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='tevox mac scripts ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.repo = self.root / 'checkout with spaces'
        for directory in ['scripts', 'src', 'docs', 'include', 'tests/data/pair']:
            shutil.copytree(REPO / directory, self.repo / directory)
        shutil.copy2(REPO / 'README.md', self.repo / 'README.md')
        shutil.copy2(REPO / 'LICENSE', self.repo / 'LICENSE')
        (self.repo / 'tests').mkdir(exist_ok=True)
        # Only check script handoff here: the real functional suite runs in Mac CI.
        (self.repo / 'tests/run_tests.sh').write_text(
            '#!/bin/sh\nset -eu\n'
            '[ "${MOCK_REGRESSION_FAIL:-0}" != 1 ] || exit 9\n'
            '"$TEVOX_BIN" --help | grep -q "^Usage:"\n'
        )
        self.tools = self.root / 'mock tools'
        self.tools.mkdir()
        for tool in ['uname', 'sw_vers', 'sysctl', 'xcrun', 'clang', 'lipo', 'codesign']:
            file = self.tools / tool
            file.write_text('#!' + sys.executable + '\n' + MOCK_TOOL)
            file.chmod(0o755)
        self.sdk = self.root / 'SDK with spaces'
        self.sdk.mkdir()
        self.log = self.root / 'invocations.jsonl'
        self.env = os.environ.copy()
        self.env.update(PATH=str(self.tools) + os.pathsep + self.env['PATH'],
                        MOCK_TOOL_LOG=str(self.log), MOCK_SDK=str(self.sdk))

    def run_script(self, script, *args, ok=True, **variables):
        env = self.env.copy()
        env.update(variables)
        result = subprocess.run(['sh', str(self.repo / 'scripts' / script), *args],
                                cwd=self.root, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if ok:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def output(self, kind='arm64'):
        return self.repo / 'build/macos' / kind / 'tevox'

    def test_native_build_quotes_paths_and_sets_sdk_and_minimum(self):
        self.run_script('build-macos.sh')
        self.assertTrue(os.access(self.output(), os.X_OK))
        calls = [json.loads(line) for line in self.log.read_text().splitlines()]
        compiles = [c for c in calls if c[0] == 'clang' and '-c' in c]
        self.assertEqual(len(compiles), len(list((self.repo / 'src').glob('*.c'))))
        for call in compiles:
            self.assertEqual(call[call.index('-arch') + 1], 'arm64')
            for flag in ['-mmacosx-version-min=11.0', '-std=c11', '-Werror', '-pedantic']:
                self.assertIn(flag, call)
            self.assertIn(str(self.sdk), call)

    def test_universal_build_contains_both_architectures(self):
        self.run_script('build-macos.sh', '--universal', MOCK_ARCH='x86_64')
        self.assertIn('# mock_archs: x86_64 arm64', self.output('universal').read_text())
        info = self.output('universal').with_name('build-info.txt').read_text()
        self.assertIn('native x86_64 slice', info)
        self.assertIn('structure verified', info)

    def test_rejects_wrong_os_old_os_and_rosetta_before_build(self):
        for variables, message in [({'MOCK_OS': 'Linux'}, 'requires macOS'),
                                   ({'MOCK_OS_VERSION': '10.15.7'}, '11.0'),
                                   ({'MOCK_ROSETTA': '1'}, 'Rosetta')]:
            with self.subTest(variables=variables):
                result = self.run_script('build-macos.sh', ok=False, **variables)
                self.assertIn(message, result.stdout)
                self.assertFalse(self.output().exists())

    def test_missing_sdk_does_not_publish_binary(self):
        self.run_script('build-macos.sh', ok=False, MOCK_SDK=str(self.root / 'missing SDK'))
        self.assertFalse(self.output().exists())

    def test_failed_build_signing_or_test_preserves_previous_binary(self):
        output = self.output()
        output.parent.mkdir(parents=True)
        output.write_bytes(b'previous known good binary')
        for variable in ['MOCK_COMPILE_FAIL', 'MOCK_SIGN_FAIL', 'MOCK_REGRESSION_FAIL']:
            with self.subTest(variable=variable):
                self.run_script('build-macos.sh', ok=False, **{variable: '1'})
                self.assertEqual(output.read_bytes(), b'previous known good binary')
                self.assertEqual(list(output.parent.parent.glob('.build.*')), [])

    def test_directory_cannot_be_used_as_binary_destination(self):
        self.output().mkdir(parents=True)
        self.run_script('build-macos.sh', ok=False)
        self.assertEqual(list(self.output().iterdir()), [])

    def test_package_has_portable_layout_and_correct_checksum(self):
        destination = self.root / 'packages with spaces'
        self.run_script('package-macos.sh', '--universal', '--revision', 'test-1',
                        '--output-dir', str(destination))
        name = 'tevox-test-1-macos-universal'
        archive = destination / (name + '.tar.gz')
        checksum = archive.with_name(archive.name + '.sha256').read_text().split()[0]
        self.assertEqual(checksum, hashlib.sha256(archive.read_bytes()).hexdigest())
        with tarfile.open(archive) as package:
            names = set(package.getnames())
            for item in ['bin/tevox', 'INSTALL.md', 'docs/build-info.txt',
                         'examples/pair/A.fa']:
                self.assertIn(name + '/' + item, names)
            self.assertTrue(package.getmember(name + '/bin/tevox').mode & 0o111)

    def test_invalid_package_label_and_failed_signing_publish_no_archive(self):
        destination = self.root / 'packages'
        self.run_script('package-macos.sh', '--revision', '../../bad', ok=False)
        self.run_script('package-macos.sh', '--output-dir', str(destination),
                        ok=False, MOCK_SIGN_FAIL='1')
        self.assertFalse(any(destination.glob('*.tar.gz')))

    def test_install_into_prefix_with_spaces_and_keep_existing_on_failure(self):
        prefix = self.root / 'user prefix'
        self.run_script('install-macos.sh', '--universal', '--prefix', str(prefix))
        binary = prefix / 'bin/tevox'
        self.assertTrue(os.access(binary, os.X_OK))
        previous = binary.read_bytes()
        self.run_script('install-macos.sh', '--prefix', str(prefix),
                        ok=False, MOCK_REGRESSION_FAIL='1')
        self.assertEqual(binary.read_bytes(), previous)
        self.assertEqual(list(binary.parent.glob('.tevox-install.*')), [])

    def test_installer_rejects_directory_and_missing_prefix(self):
        prefix = self.root / 'user prefix'
        (prefix / 'bin/tevox').mkdir(parents=True)
        self.run_script('install-macos.sh', '--prefix', str(prefix), ok=False)
        self.assertEqual(list((prefix / 'bin/tevox').iterdir()), [])
        self.run_script('install-macos.sh', '--prefix', ok=False)


if __name__ == '__main__':
    print('Testing simulated Mac tool orchestration; this is NOT native macOS validation.', flush=True)
    unittest.main(verbosity=2)

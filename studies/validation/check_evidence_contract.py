#!/usr/bin/env python3
"""Stream per-row evidence invariants without loading a whole genome graph.

This supplements counts/hash checks; it does not replace the complete relational
schema validator used on bounded fixtures and the Arabidopsis pair.
"""
import argparse
import collections
import csv
import gzip
import json
import math
from pathlib import Path


def check(prefix):
    prefix = Path(prefix)
    run = json.loads(Path(str(prefix) + '.run.json').read_text())
    path = Path(str(prefix) + '.evidence.tsv.gz')
    if not path.exists():
        path = Path(str(prefix) + '.evidence.tsv')
    opener = gzip.open if path.suffix == '.gz' else open
    methods = collections.Counter()
    count = 0
    with opener(path, 'rt') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        ix = {name: header.index(name) for name in [
            'schema_version', 'local_identity', 'identity_method',
            'alignment_identity', 'alignment_identity_method', 'mapq',
            'mapq_status', 'legacy_state', 'mapping_confidence', 'origin',
            'dependency', 'nearby_candidate_count', 'retained_candidate_count',
            'graph_candidate_count', 'provider']}
        for count, row in enumerate(reader, 1):
            def require(condition, message):
                if not condition:
                    raise ValueError('{}:{}: {}'.format(path, count + 1, message))
            require(len(row) == len(header), 'column count')
            get = lambda name: row[ix[name]]
            require(get('schema_version') == run['schema_version'], 'schema version')
            for value, method in [('local_identity', 'identity_method'),
                                  ('alignment_identity', 'alignment_identity_method')]:
                missing = get(value) == '.'
                require(missing == (get(method) == 'MISSING'), value + ' availability')
                if not missing:
                    number = float(get(value))
                    require(math.isfinite(number) and 0 <= number <= 1, value + ' range')
            methods[get('identity_method')] += 1
            require(get('mapq_status') in {'MISSING_255', 'NOT_PROVIDED', 'OBSERVED'}, 'mapq status')
            require((get('mapq') == '.') == (get('mapq_status') != 'OBSERVED'), 'mapq availability')
            require((get('origin'), get('dependency')) in {
                ('DERIVED_REVERSE', 'DERIVED_SAME_GROUP'), ('NATIVE', 'INDEPENDENT')}, 'dependency')
            if get('legacy_state') == 'EMPTY_SITE_CONFIRMED':
                require(get('provider') == 'PAF' and get('local_identity') != '.' and
                        get('mapping_confidence') != 'NOT_ESTABLISHED', 'empty-site evidence')
            nearby = int(get('nearby_candidate_count'))
            require(nearby >= int(get('retained_candidate_count')) >= 0 and
                    nearby >= int(get('graph_candidate_count')) >= 0, 'candidate counts')
    if count != run['counts']['evidence_observations']:
        raise ValueError('evidence count differs from run manifest')
    return {'status': 'PASS', 'path': str(path), 'rows': count,
            'local_identity_methods': dict(methods),
            'scope': 'streamed per-row evidence invariants; not full relational schema'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefixes', nargs='+', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    results = []
    for prefix in args.prefixes:
        result = check(prefix)
        results.append(result)
        print(prefix, result['status'], result['rows'], flush=True)
    args.output.write_text(json.dumps({'status': 'PASS', 'runs': results}, indent=2) + '\n')

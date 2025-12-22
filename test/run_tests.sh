#!/bin/bash

# Test script for TE Comparator

echo "=== TE Comparator Test Script ==="
echo

# Check if the executable exists
if [ ! -f "./tevox" ]; then
    echo "Error: tevo executable not found. Please run 'make' first."
    exit 1
fi

# Test 1: Basic functionality test
echo "Test 1: Basic functionality test"
echo "Running: ./tevox test_data/synteny_example.txt test_data/genome1_te.gff3 test_data/genome2_te.bed -v"
echo

./tevox test_data/synteny_example.txt test_data/genome1_te.gff3 test_data/genome2_te.bed -v

if [ $? -eq 0 ]; then
    echo "✓ Test 1 passed"
else
    echo "✗ Test 1 failed"
fi

echo
echo "====================================="
echo

# Test 2: Test with custom output prefix
echo "Test 2: Custom output prefix test"
echo "Running: ./tevox test_data/synteny_example.txt test_data/genome1_te.gff3 test_data/genome2_te.bed -o test_output"
echo

./tevox test_data/synteny_example.txt test_data/genome1_te.gff3 test_data/genome2_te.bed -o test_output

if [ $? -eq 0 ]; then
    echo "✓ Test 2 passed"
    
    # Check if output files were created
    if [ -f "test_output_genome1_unique.txt" ] && [ -f "test_output_genome2_unique.txt" ]; then
        echo "✓ Output files created successfully"
        echo "  - test_output_genome1_unique.txt"
        echo "  - test_output_genome2_unique.txt"
    else
        echo "✗ Output files not found"
    fi
else
    echo "✗ Test 2 failed"
fi

echo
echo "====================================="
echo

# Test 3: Help message test
echo "Test 3: Help message test"
echo "Running: ./tevox --help"
echo

./tevox --help

if [ $? -eq 0 ]; then
    echo "✓ Test 3 passed"
else
    echo "✗ Test 3 failed"
fi

echo
echo "====================================="
echo

# Test 4: Error handling test (invalid file)
echo "Test 4: Error handling test"
echo "Running: ./tevox nonexistent.txt test_data/genome1_te.gff3 test_data/genome2_te.bed"
echo

./tevox nonexistent.txt test_data/genome1_te.gff3 test_data/genome2_te.bed

if [ $? -ne 0 ]; then
    echo "✓ Test 4 passed (correctly handled error)"
else
    echo "✗ Test 4 failed (should have failed)"
fi

echo
echo "=== Test Summary ==="
echo "All tests completed. Check the output above for any failures."
echo "Generated files can be found in the current directory."
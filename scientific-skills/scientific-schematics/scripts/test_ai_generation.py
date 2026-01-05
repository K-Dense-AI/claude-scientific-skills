#!/usr/bin/env python3
"""
Test script to verify AI generation implementation.

This script performs dry-run tests without making actual API calls.
It verifies:
1. Script structure and imports
2. Class initialization
3. Method signatures
4. Error handling
5. Command-line interface

Usage:
    python test_ai_generation.py
"""

import sys
import os
from pathlib import Path

# Add scripts directory to path
scripts_dir = Path(__file__).parent / "scripts"
sys.path.insert(0, str(scripts_dir))

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    try:
        from generate_schematic_ai import ScientificSchematicGenerator
        print("✓ generate_schematic_ai imports successfully")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False

def test_class_structure():
    """Test class initialization and structure."""
    print("\nTesting class structure...")
    try:
        from generate_schematic_ai import ScientificSchematicGenerator
        
        # Test initialization with dummy key
        generator = ScientificSchematicGenerator(api_key="test_key", verbose=False)
        print("✓ Class initializes successfully")
        
        # Check required methods exist
        required_methods = [
            'generate_image',
            'review_image',
            'improve_prompt',
            'generate_iterative'
        ]
        
        for method in required_methods:
            if not hasattr(generator, method):
                print(f"✗ Missing method: {method}")
                return False
            print(f"✓ Method exists: {method}")
        
        # Check attributes
        if not hasattr(generator, 'api_key'):
            print("✗ Missing attribute: api_key")
            return False
        print("✓ Attribute exists: api_key")
        
        if not hasattr(generator, 'image_model'):
            print("✗ Missing attribute: image_model")
            return False
        print(f"✓ Image model: {generator.image_model}")
        
        if not hasattr(generator, 'review_model'):
            print("✗ Missing attribute: review_model")
            return False
        print(f"✓ Review model: {generator.review_model}")
        
        return True
    except Exception as e:
        print(f"✗ Class structure test failed: {e}")
        return False

def test_error_handling():
    """Test error handling for missing API key."""
    print("\nTesting error handling...")
    try:
        from generate_schematic_ai import ScientificSchematicGenerator
        
        # Clear environment variable
        old_key = os.environ.get("OPENROUTER_API_KEY")
        if old_key:
            del os.environ["OPENROUTER_API_KEY"]
        
        # Try to initialize without key
        try:
            generator = ScientificSchematicGenerator()
            print("✗ Should have raised ValueError for missing API key")
            return False
        except ValueError as e:
            if "OPENROUTER_API_KEY" in str(e):
                print("✓ Correctly raises ValueError for missing API key")
            else:
                print(f"✗ Wrong error message: {e}")
                return False
        
        # Restore environment variable
        if old_key:
            os.environ["OPENROUTER_API_KEY"] = old_key
        
        return True
    except Exception as e:
        print(f"✗ Error handling test failed: {e}")
        return False

def test_wrapper_script():
    """Test wrapper script structure."""
    print("\nTesting wrapper script...")
    try:
        import generate_schematic
        print("✓ generate_schematic imports successfully")
        
        # Check main functions exist
        if not hasattr(generate_schematic, 'main'):
            print("✗ Missing function: main")
            return False
        print("✓ Function exists: main")
        
        return True
    except Exception as e:
        print(f"✗ Wrapper script test failed: {e}")
        return False

def test_prompt_engineering():
    """Test prompt construction."""
    print("\nTesting prompt engineering...")
    try:
        from generate_schematic_ai import ScientificSchematicGenerator
        
        generator = ScientificSchematicGenerator(api_key="test_key", verbose=False)
        
        # Test improve_prompt method
        original = "Create a flowchart"
        critique = "Add more spacing between boxes"
        improved = generator.improve_prompt(original, critique, 2)
        
        if not improved:
            print("✗ improve_prompt returned empty string")
            return False
        
        if original not in improved:
            print("✗ Improved prompt doesn't include original")
            return False
        
        if critique not in improved:
            print("✗ Improved prompt doesn't include critique")
            return False
        
        if "ITERATION 2" not in improved:
            print("✗ Improved prompt doesn't include iteration number")
            return False
        
        print("✓ Prompt engineering works correctly")
        print(f"  Original length: {len(original)} chars")
        print(f"  Improved length: {len(improved)} chars")
        
        return True
    except Exception as e:
        print(f"✗ Prompt engineering test failed: {e}")
        return False

def test_file_paths():
    """Test that all required files exist."""
    print("\nTesting file structure...")
    
    base_dir = Path(__file__).parent
    required_files = [
        "scripts/generate_schematic_ai.py",
        "scripts/generate_schematic.py",
        "SKILL.md",
        "README.md"
    ]
    
    all_exist = True
    for file_path in required_files:
        full_path = base_dir / file_path
        if full_path.exists():
            print(f"✓ {file_path}")
        else:
            print(f"✗ Missing: {file_path}")
            all_exist = False
    
    return all_exist

def main():
    """Run all tests."""
    print("="*60)
    print("Scientific Schematics AI Generation - Verification Tests")
    print("="*60)
    
    tests = [
        ("File Structure", test_file_paths),
        ("Imports", test_imports),
        ("Class Structure", test_class_structure),
        ("Error Handling", test_error_handling),
        ("Wrapper Script", test_wrapper_script),
        ("Prompt Engineering", test_prompt_engineering),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n✗ Test '{test_name}' crashed: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n✓ All tests passed! Implementation verified.")
        print("\nNext steps:")
        print("1. Set OPENROUTER_API_KEY environment variable")
        print("2. Test with actual API call:")
        print("   python scripts/generate_schematic.py 'test diagram' -o test.png")
        return 0
    else:
        print(f"\n✗ {total - passed} test(s) failed. Please review errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())


#!/usr/bin/env python3
"""
validate_filter.py - Diagnostic script for mzML Filter workspace environment validation.
Checks imports, module availability, and file access to ensure all components are properly installed.
Provides detailed status reporting for PyOpenMS, Streamlit, and custom filtering modules.
"""

import sys
from pathlib import Path

def check_imports():
    """Check that all required modules can be imported."""
    print("🔍 Checking imports...")
    
    try:
        import pyopenms
        print("  ✅ PyOpenMS")
    except ImportError as e:
        print(f"  ❌ PyOpenMS: {e}")
        return False
    
    try:
        import streamlit
        print("  ✅ Streamlit")
    except ImportError as e:
        print(f"  ❌ Streamlit: {e}")
        return False
    
    try:
        import logging
        print("  ✅ logging")
    except ImportError as e:
        print(f"  ❌ logging: {e}")
        return False
    
    return True


def check_module_import():
    """Check that mzml_processor module can be imported."""
    print("\n🔍 Checking mzml_processor module...")
    
    try:
        core_dir = Path(__file__).parent.parent / "core"
        sys.path.insert(0, str(core_dir))
        
        from mzml_processor import (
            extract_peptide_relevant_spectra,
            filter_mzml,
            process_idxml_batch
        )
        print("  ✅ All functions imported successfully")
        return True
    except Exception as e:
        print(f"  ❌ Import failed: {e}")
        return False


def check_files():
    """Check that all required files exist."""
    print("\n🔍 Checking files...")
    
    base_dir = Path(__file__).parent.parent
    
    required_files = [
        "core/mzml_processor.py",
        "cli/filter_mzml.py",
        "cli/run_mzml_filter.sh",
        "workspaces/dataset_filter.py",
        "docs/MZML_FILTER.md",
        "docs/IMPLEMENTATION_SUMMARY.md",
        "examples/mzml_filter_examples.py"
    ]
    
    all_exist = True
    for file_path in required_files:
        full_path = base_dir / file_path
        if full_path.exists():
            print(f"  ✅ {file_path}")
        else:
            print(f"  ❌ {file_path} (missing)")
            all_exist = False
    
    return all_exist


def check_executables():
    """Check that scripts are executable."""
    print("\n🔍 Checking executables...")
    
    base_dir = Path(__file__).parent.parent
    
    scripts = [
        "bin/run_mzml_filter.sh",
        "bin/mzml_filter_cli.py"
    ]
    
    all_executable = True
    for script_path in scripts:
        full_path = base_dir / script_path
        if full_path.exists():
            import os
            if os.access(full_path, os.X_OK):
                print(f"  ✅ {script_path}")
            else:
                print(f"  ⚠️  {script_path} (not executable, run: chmod +x {full_path})")
                all_executable = False
        else:
            print(f"  ❌ {script_path} (missing)")
            all_executable = False
    
    return all_executable


def main():
    print("="*70)
    print("mzML Filter - Validation Script")
    print("="*70)
    
    results = {
        "Imports": check_imports(),
        "Module": check_module_import(),
        "Files": check_files(),
        "Executables": check_executables()
    }
    
    print("\n" + "="*70)
    print("Validation Results")
    print("="*70)
    
    for check, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{check:20s} {status}")
    
    print("="*70)
    
    if all(results.values()):
        print("\n🎉 All checks passed! mzML Filter is ready to use.\n")
        print("Quick start:")
        print("  Web Interface: ./bin/run_mzml_filter.sh")
        print("  CLI:           ./bin/mzml_filter_cli.py --help")
        print("  Examples:      python examples/mzml_filter_examples.py")
        return 0
    else:
        print("\n⚠️  Some checks failed. Please review the errors above.\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
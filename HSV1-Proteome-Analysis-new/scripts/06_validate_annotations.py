"""
Script 06: Validate Annotations (Wrapper for Biological Temporal Annotation Audit)

Executes the biological temporal annotation audit to verify all 74 unique proteins
against peer-reviewed experimental virology literature.
"""

import os
import sys
import importlib.util

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    audit_script = os.path.join(script_dir, "06_audit_annotations.py")
    
    spec = importlib.util.spec_from_file_location("audit_module", audit_script)
    audit_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit_module)
    
    audit_module.run_audit()

if __name__ == "__main__":
    main()

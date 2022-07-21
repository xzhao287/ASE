"""Check that plain cli doesn't execute too many imports."""
from ase.utils.checkimports import check_imports


def test_imports():
    forbidden_modules = [
        'gpaw',  # external
        'scipy',  # large
        'ase.io.formats',  # possibly slow external formats
        'ase.calculators.(?!names).*',  # any calculator
    ]
    check_imports("from ase.cli.main import main; main(args=[])",
                  forbidden_modules=forbidden_modules,
                  max_module_count=270)

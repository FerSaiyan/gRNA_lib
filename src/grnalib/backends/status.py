from __future__ import annotations

from importlib.util import find_spec
from shutil import which


def backend_status() -> dict:
    """Return lightweight capability information without importing heavy ML stacks."""
    return {
        'rs3': {
            'available': find_spec('rs3') is not None,
            'kind': 'python',
            'purpose': 'SpCas9 on-target activity',
        },
        'deepprime_genet': {
            'available': find_spec('genet') is not None,
            'kind': 'python',
            'purpose': 'learned prime-edit efficiency',
            'note': 'GenET currently fits best in a Python 3.10 environment because of its TensorFlow metadata.',
        },
        'crisprware': {
            'available': which('crisprware') is not None,
            'kind': 'external-command',
            'purpose': 'indexed off-target specificity using crispr-ots or GuideScan2 index',
            'license_note': (
                'CRISPRware is externally installed and is not vendored by gRNA Library. '
                'Its repository currently uses the UC Santa Cruz Noncommercial License; '
                'commercial users should review or obtain the appropriate license.'
            ),
        },
    }

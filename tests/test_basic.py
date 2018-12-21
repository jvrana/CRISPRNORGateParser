import pytest
from norgateparser import parse_nor_gate_name_to_sequence



@pytest.mark.parametrize('name', [
    'pMOD-LTR1-Nat-pGRR-W5W8-URGR-F1',
    'pMOD6-pGRR-F1-yeGFP',
    'PMOD6-PGRR-F1-yeGFP',
    'pMOD8-pGRR-W5W8-iRGR-W36',
    'pMOD8-pGRR-W8-iRGR-W36',
])
def test_basic(name):
    seq = parse_nor_gate_name_to_sequence(name)
    assert seq is not None
    seq.print()

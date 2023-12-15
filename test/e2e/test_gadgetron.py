
import pytest
from pathlib import Path

def cases():
    cases_path = Path(__file__).parent.parent.joinpath("integration/cases")
    return cases_path.glob("*.cfg")    

@pytest.mark.parametrize("test_case", cases())
def test_gadgetron(test_case, test_orchestrator_generator):
    orchestrator = test_orchestrator_generator(test_case)
    orchestrator.init()
    orchestrator.validate_requirements()
    orchestrator.run()
    orchestrator.end()

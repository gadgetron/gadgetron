
import pytest
from typing import List

class TestStep:
    def init(self):
        pass

    def run(self):
        pass

    def validate(self):
        pass

    def end(self):
        pass

class TestOrchestrator:
    def __init__(self) -> None:
        self.steps: List[TestStep] = []

    def init(self):
        pass

    def validate_requirements(self):
        pass

    def run(self):
        for step in self.steps:
            step.init()
            step.run()
            step.validate()
            step.end()

    def end(self):
        pass


@pytest.fixture(scope="module")
def test_orchestrator_generator():

    def generate_test_orchestrator(config_file: str):
        return TestOrchestrator()
    
    return generate_test_orchestrator


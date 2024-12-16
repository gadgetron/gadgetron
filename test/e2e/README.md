## End-to-End Tests

The end-to-end tests use `Pytest` for test management and invocation.

When a test is run for the first time, the corresponding test data will be downloaded *automatically*.

Examples:

- `pytest`: Run all compatible end-to-end tests
- `pytest -rA -s`: Print test steps as they are run
- `pytest -k simple`: Run all tests with "simple" in the name
- `pytest --tags fast`: Run all tests with the tag "fast"
- `pytest --ignore-requirements=gpu_support`: Run tests requiring GPUs, even if GPU is not available
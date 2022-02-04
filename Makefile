PROJECT_NAME = "diascope"

.PHONY: venv test

venv:
	poetry install

kernel: venv
	poetry run python -m ipykernel install --user --name ${PROJECT_NAME}

test: venv
	poetry run python -m pytest --durations=0 -s $(FILTER)
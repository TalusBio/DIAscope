.PHONY: venv test

venv:
	poetry install

test: venv
	poetry run python -m pytest --durations=0 -s $(FILTER)
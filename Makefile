# Offline targets: none of these need the data, the containers, or a network connection.
.PHONY: test numbers guide

test:            ## check README/docs against results/ and the code
	python3 -m pytest -q tests

numbers:         ## print the headline numbers straight from results/
	python3 scripts/headline_numbers.py

guide:           ## rebuild docs/interview-guide.pdf from results/ (needs reportlab)
	python3 scripts/build_interview_guide.py

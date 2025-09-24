# Makefile for JointODE presentation
.PHONY: all clean presenter standard print handout

# Default - build all versions
all: standard presenter print handout clean

# Standard presentation (no notes)
standard:
	@echo "Building standard presentation..."
	@pdflatex -interaction=nonstopmode presentation.tex > /dev/null 2>&1 || true
	@pdflatex -interaction=nonstopmode presentation.tex > /dev/null 2>&1 || true
	@echo "✓ Created presentation.pdf"

# Presenter mode (with notes on second screen)
presenter:
	@echo "Building presenter version..."
	@cp presentation.tex temp.tex
	@sed -i '' 's|%\\setbeameroption{show notes on second screen=right}|\\setbeameroption{show notes on second screen=right}|' temp.tex
	@pdflatex -interaction=nonstopmode temp.tex > /dev/null 2>&1 || true
	@pdflatex -interaction=nonstopmode temp.tex > /dev/null 2>&1 || true
	@mv temp.pdf presentation_presenter.pdf
	@rm -f temp.*
	@echo "✓ Created presentation_presenter.pdf"

# Print version (slides + notes)
print:
	@echo "Building print version..."
	@cp presentation.tex temp.tex
	@sed -i '' 's|%\\setbeameroption{show notes}|\\setbeameroption{show notes}|' temp.tex
	@pdflatex -interaction=nonstopmode temp.tex > /dev/null 2>&1 || true
	@pdflatex -interaction=nonstopmode temp.tex > /dev/null 2>&1 || true
	@mv temp.pdf presentation_print.pdf
	@rm -f temp.*
	@echo "✓ Created presentation_print.pdf"

# Handout (4 slides per page)
handout:
	@echo "Building handout..."
	@cp presentation.tex temp.tex
	@sed -i '' 's|%\\pgfpagesuselayout{4 on 1}|\\pgfpagesuselayout{4 on 1}|' temp.tex
	@pdflatex -interaction=nonstopmode temp.tex > /dev/null 2>&1 || true
	@pdflatex -interaction=nonstopmode temp.tex > /dev/null 2>&1 || true
	@mv temp.pdf presentation_handout.pdf
	@rm -f temp.*
	@echo "✓ Created presentation_handout.pdf"

# Clean auxiliary files
clean:
	@rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *.synctex.gz
	@echo "✓ Cleaned"
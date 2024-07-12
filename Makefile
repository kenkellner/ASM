SRC=src
BUILD=build

RMD_IN = $(wildcard $(SRC)/*.Rmd)
RMD_OUT := $(patsubst $(SRC)/%.Rmd,$(BUILD)/%.html,$(RMD_IN))
CODE_IN = $(wildcard $(SRC)/code/*.R)

all: $(RMD_OUT)
	@mkdir -p build
	@echo "Building index"
	@Rscript build_Rmd.R $(SRC)/index.Rmd $(BUILD)/index.html > /dev/null 2>&1
	@echo "Done"

code: $(CODE_IN) $(SRC)/code/code_page_template.Rmd
	@mkdir -p $(BUILD)/code
	@cp $(SRC)/code/*.html $(BUILD)/code/

deploy:
	@rsync -r --progress --delete --update build/ \
		lemur:/var/www/kenkellner.com/ASM/

$(BUILD)/%.html: $(SRC)/%.Rmd $(SRC)/_navbar.yml $(SRC)/style.css 
	Rscript build_Rmd.R $< $@

clean:
	rm -f build/*

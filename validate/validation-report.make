INPUT = $(wildcard identicals.??)
OUTPUT = $(INPUT:identicals.%=counts.%)

counts.%: identicals.%
	./validate-all.sh $^ > $@

validation-report.tsv: $(OUTPUT)
	cat $^ > $@

all: validation-report.tsv


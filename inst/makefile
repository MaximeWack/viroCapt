DATADIR=fastqs

QCDIR=01_QC
GLOBALDIR=02_Global
LOCALDIR=03_Local
FASTADIR=04_Fasta
BLATDIR=05_Blat
VISUDIR=06_Visu

THREADS=2
BT_OPTIONS_GLOBAL=--no-unal --very-sensitive
BT_OPTIONS_LOCAL=--local --no-unal --very-sensitive-local
BLAT_OPTIONS=-minScore=25 -minIdentity=90 -noHead

LOCAL_PLOT_LIMIT=5
FINAL_PLOT_LIMIT=1

red="\033[0;31m"
end="\033[0m"

fastqs=$(shell \
	find $(DATADIR) -maxdepth 1 -type f -name '*.fastq*' \
	| sed -r 's%\.(R[[:digit:]]|U)\.fastq.*%%' \
	| sed -r 's%$(DATADIR)\/%%' \
	| sort \
	| uniq)

All: QC Local_plot Consensus Final_plot Visu

# Create an R object to use in the visualisation
Visu: $(patsubst %, $(DATADIR)/$(VISUDIR)/%.rds, $(fastqs))

# Plot the final Plots
Final_plot: $(patsubst %, $(DATADIR)/$(BLATDIR)/%.png, $(fastqs))

# Summarise the cleaned blat file
Summarise_blat: $(patsubst %, $(DATADIR)/$(BLATDIR)/%_summary.tsv, $(fastqs))

# Clean the blat file
Clean_blat: $(patsubst %, $(DATADIR)/$(BLATDIR)/%_cleaned.tsv, $(fastqs))

# Blat the fasta files
Blat: $(patsubst %, $(DATADIR)/$(BLATDIR)/%.tsv, $(fastqs))

# Extract fasta files for non-viral sequences
Extract_fasta: $(patsubst %, $(DATADIR)/$(FASTADIR)/%.fa, $(fastqs))

# Generate consensus
Consensus: $(patsubst %, $(DATADIR)/$(LOCALDIR)/%_consensus.fa, $(fastqs))

# Plots for local alignment
Local_plot: $(patsubst %, $(DATADIR)/$(LOCALDIR)/%.png, $(fastqs))

# # Sequencing profile generation for plots
Profile: $(patsubst %, $(DATADIR)/$(LOCALDIR)/%.rds, $(fastqs))

# Local alignment
Local: $(patsubst %, $(DATADIR)/$(LOCALDIR)/%.sam, $(fastqs))

# Selection of genotypes
Selection: $(patsubst %, $(DATADIR)/$(GLOBALDIR)/%.tsv.sel, $(fastqs))

# Genotyping
Genotyping: $(patsubst %, $(DATADIR)/$(GLOBALDIR)/%.tsv, $(fastqs))

# Global alignment
# Global: $(patsubst %, $(DATADIR)/$(GLOBALDIR)/%.sam, $(fastqs))

# QC
QC: $(patsubst %, $(DATADIR)/$(QCDIR)/%.png, $(fastqs))

# QC_sam
# QC_sam: $(patsubst %, $(DATADIR)/$(QCDIR)/%.sam, $(fastqs))

# Create an R object to use in the visualisation
$(DATADIR)/$(VISUDIR)/%.rds: $(DATADIR)/$(BLATDIR)/%_summary.tsv $(DATADIR)/$(LOCALDIR)/%.rds $(DATADIR)/$(VISUDIR)
	@echo -e $(red) "Creating the visualisation object for $*" $(end)
	Rscript -e 'viroCapt::create_visu("$(DATADIR)/$(LOCALDIR)/$*", "$(DATADIR)/$(BLATDIR)/$*", "$@")'

# Plot a final plot
$(DATADIR)/$(BLATDIR)/%.png: $(DATADIR)/$(BLATDIR)/%_summary.tsv $(DATADIR)/$(LOCALDIR)/%.rds
	@echo -e $(red) "Plotting the final plot for $*" $(end)
	Rscript -e 'viroCapt::plot_final("$(DATADIR)/$(LOCALDIR)/$*", "$(DATADIR)/$(BLATDIR)/$*", limit = $(FINAL_PLOT_LIMIT))'

# Summarise a cleaned blat file
$(DATADIR)/$(BLATDIR)/%_summary.tsv: $(DATADIR)/$(BLATDIR)/%_cleaned.tsv
	@echo -e $(red) "Summarizing cleaned blat file for $*" $(end)
	Rscript -e 'viroCapt::summarised_blat("$<", "$@")'

# Clean a blat file
$(DATADIR)/$(BLATDIR)/%_cleaned.tsv: $(DATADIR)/$(BLATDIR)/%.tsv $(DATADIR)/$(FASTADIR)/%.fa
	@echo -e $(red) "Cleaning blat file for $*" $(end)
	Rscript -e 'viroCapt::cleaned_blat("$(word 1, $^)", "$(word 2, $^)", "$@")'

# Blat the fasta files
$(DATADIR)/$(BLATDIR)/%.tsv: $(DATADIR)/$(FASTADIR)/%.fa | $(DATADIR)/$(BLATDIR) hg19.2bit
	@echo -e $(red) "Blat non-viral sequences for $*" $(end)
	blat hg19.2bit $< $@ $(BLAT_OPTIONS)

# Extract non-viral sequences as fasta files
$(DATADIR)/$(FASTADIR)/%.fa: $(DATADIR)/$(LOCALDIR)/%.sam | $(DATADIR)/$(FASTADIR)
	@echo -e $(red) "Extracting non-viral sequences for $*" $(end)
	Rscript -e 'viroCapt::extract_fasta("$<", "$@", $(THREADS))'

# Consensus files
$(DATADIR)/$(LOCALDIR)/%_consensus.fa: $(DATADIR)/$(LOCALDIR)/%.rds
	@echo -e $(red) "Writing computed consensus for $*" $(end)
	Rscript -e 'viroCapt::write_consensus("$(DATADIR)/$(LOCALDIR)/$*")'

# Local alignment plot
$(DATADIR)/$(LOCALDIR)/%.png: $(DATADIR)/$(LOCALDIR)/%.rds
	@echo -e $(red) "Plotting local alignment for $*" $(end)
	Rscript -e 'viroCapt::plot_depth("$(DATADIR)/$(LOCALDIR)/$*", limit = $(LOCAL_PLOT_LIMIT))'

# Profile
$(DATADIR)/$(LOCALDIR)/%.rds: $(DATADIR)/$(LOCALDIR)/%.sam
	@echo -e $(red) "Computing sequencing profile and consensus sequences for $*" $(end)
	Rscript -e 'viroCapt::create_profile("$(DATADIR)/$(LOCALDIR)/$*", $(THREADS), consensus = T)'

## Local alignment pairend gzipped
$(DATADIR)/$(LOCALDIR)/%.sam: $(DATADIR)/$(GLOBALDIR)/%.tsv.sel $(DATADIR)/%.R1.fastq.gz $(DATADIR)/%.R2.fastq.gz | $(DATADIR)/$(LOCALDIR) Refs
	for genotype in `cat $(DATADIR)/$(GLOBALDIR)/$*.tsv.sel`; do \
		echo -e $(red) "Local alignment on $$genotype for $*" $(end) ; \
		bowtie2 -x refs/$$genotype -1 $(DATADIR)/$*.R1.fastq.gz -2 $(DATADIR)/$*.R2.fastq.gz $(BT_OPTIONS_LOCAL) -p $(THREADS) \
			| samtools sort -@ $(THREADS) -O sam - \
			| samtools rmdup - - --output-fmt sam \
			| samtools view >> $@ ; \
	done

## Local alignment pairend
$(DATADIR)/$(LOCALDIR)/%.sam: $(DATADIR)/$(GLOBALDIR)/%.tsv.sel $(DATADIR)/%.R1.fastq $(DATADIR)/%.R2.fastq | $(DATADIR)/$(LOCALDIR) Refs
	for genotype in `cat $(DATADIR)/$(GLOBALDIR)/$*.tsv.sel`; do \
		echo -e $(red) "Local alignment on $$genotype for $*" $(end) ; \
		bowtie2 -x refs/$$genotype -1 $(DATADIR)/$*.R1.fastq -2 $(DATADIR)/$*.R2.fastq $(BT_OPTIONS_LOCAL) -p $(THREADS) \
			| samtools sort -@ $(THREADS) -O sam - \
			| samtools rmdup - - --output-fmt sam \
			| samtools view >> $@ ; \
	done

## Local alignment single gzipped
$(DATADIR)/$(LOCALDIR)/%.sam: $(DATADIR)/$(GLOBALDIR)/%.tsv.sel $(DATADIR)/%.U.fastq.gz | $(DATADIR)/$(LOCALDIR) Refs
	for genotype in `cat $(DATADIR)/$(GLOBALDIR)/$*.tsv.sel`; do \
		echo -e $(red) "Local alignment on $$genotype for $*" $(end) ; \
		bowtie2 -x refs/$$genotype -U $(DATADIR)/$*.U.fastq.gz $(BT_OPTIONS_LOCAL) -p $(THREADS) \
			| samtools sort -@ $(THREADS) -O sam - \
			| samtools view >> $@ ; \
	done

## Local alignment single
$(DATADIR)/$(LOCALDIR)/%.sam: $(DATADIR)/$(GLOBALDIR)/%.tsv.sel $(DATADIR)/%.U.fastq | $(DATADIR)/$(LOCALDIR) Refs
	for genotype in `cat $(DATADIR)/$(GLOBALDIR)/$*.tsv.sel`; do \
		echo -e $(red) "Local alignment on $$genotype for $*" $(end) ; \
		bowtie2 -x refs/$$genotype -U $(DATADIR)/$*.U.fastq $(BT_OPTIONS_LOCAL) -p $(THREADS) \
			| samtools sort -@ $(THREADS) -O sam - \
			| samtools view >> $@ ; \
	done

## Genotyping

# Select best variant for each genotype
$(DATADIR)/$(GLOBALDIR)/%.tsv.sel: $(DATADIR)/$(GLOBALDIR)/%.tsv
	@echo -e $(red) "Selecting genotype+variant for $*" $(end)
	sed -r 's/^ +//' $< \
		| sed -r 's/ /\t/' \
		| sed -r 's/_/\t/' \
		| sort -k2,2 -k1,1nr \
		| sort -k2,2 -u \
		| sort -nr  \
		| sed -r 's/^.+\tH/H/' \
		| sed -r 's/\t/_/' > $@

# Extract genotypes from global alignment
$(DATADIR)/$(GLOBALDIR)/%.tsv: $(DATADIR)/$(GLOBALDIR)/%.sam
	@echo -e $(red) "Genotyping $*" $(end)
	cut $< -f 3 \
		| uniq -c \
		| sort -nr > $@

## Global alignment

# Global alignment pairend gzipped
$(DATADIR)/$(GLOBALDIR)/%.sam: $(DATADIR)/%.R1.fastq.gz $(DATADIR)/%.R2.fastq.gz | $(DATADIR)/$(GLOBALDIR) Refs
	@echo -e $(red) "Running pair-end global alignment for $*" $(end)
	bowtie2 -x refs/HPVs -1 $(DATADIR)/$*.R1.fastq.gz -2 $(DATADIR)/$*.R2.fastq.gz $(BT_OPTIONS_GLOBAL) -p $(THREADS) \
		| samtools sort -@ $(THREADS) -O sam - \
		| samtools view > $@

# Global alignment pairend
$(DATADIR)/$(GLOBALDIR)/%.sam: $(DATADIR)/%.R1.fastq $(DATADIR)/%.R2.fastq | $(DATADIR)/$(GLOBALDIR) Refs
	@echo -e $(red) "Running pair-end global alignment for $*" $(end)
	bowtie2 -x refs/HPVs -1 $(DATADIR)/$*.R1.fastq -2 $(DATADIR)/$*.R2.fastq $(BT_OPTIONS_GLOBAL) -p $(THREADS) \
		| samtools sort -@ $(THREADS) -O sam - \
		| samtools view > $@

# Global alignment single gzipped
$(DATADIR)/$(GLOBALDIR)/%.sam: $(DATADIR)/%.U.fastq.gz | $(DATADIR)/$(GLOBALDIR) Refs
	@echo -e $(red) "Running single read global alignment for $*" $(end)
	bowtie2 -x refs/HPVs -U $(DATADIR)/$*.U.fastq.gz $(BT_OPTIONS_GLOBAL) -p $(THREADS) \
		| samtools sort -@ $(THREADS) -O sam - \
		| samtools view > $@

# Global alignment single
$(DATADIR)/$(GLOBALDIR)/%.sam: $(DATADIR)/%.U.fastq | $(DATADIR)/$(GLOBALDIR) Refs
	@echo -e $(red) "Running single read global alignment for $*" $(end)
	bowtie2 -x refs/HPVs -U $(DATADIR)/$*.U.fastq $(BT_OPTIONS_GLOBAL) -p $(THREADS) \
		| samtools sort -@ $(THREADS) -O sam - \
		| samtools view > $@

## QC

# Quality control plot
$(DATADIR)/$(QCDIR)/%.png: $(DATADIR)/$(QCDIR)/%.rds
	@echo -e $(red) "Plotting quality control for $*" $(end)
	Rscript -e 'viroCapt::plot_depth("$(DATADIR)/$(QCDIR)/$*")'

# Quality control profile
$(DATADIR)/$(QCDIR)/%.rds: $(DATADIR)/$(QCDIR)/%.sam
	@echo -e $(red) "Computing sequencing profile for $*" $(end)
	Rscript -e 'viroCapt::create_profile("$(DATADIR)/$(QCDIR)/$*", $(THREADS))'

# Quality control pairend gzipped
$(DATADIR)/$(QCDIR)/%.sam: $(DATADIR)/%.R1.fastq.gz $(DATADIR)/%.R2.fastq.gz | $(DATADIR)/$(QCDIR) Refs
	@echo -e $(red) "Running pair-end quality control for $*" $(end)
	bowtie2 -x refs/ctrl_capture -1 $(DATADIR)/$*.R1.fastq.gz -2 $(DATADIR)/$*.R2.fastq.gz $(BT_OPTIONS_GLOBAL) -p $(THREADS) \
		| samtools sort -@ $(THREADS) -O sam - \
		| samtools view > $@

# Quality control pairend
$(DATADIR)/$(QCDIR)/%.sam: $(DATADIR)/%.R1.fastq $(DATADIR)/%.R2.fastq | $(DATADIR)/$(QCDIR) Refs
	@echo -e $(red) "Running pair-end quality control for $*" $(end)
	bowtie2 -x refs/ctrl_capture -1 $(DATADIR)/$*.R1.fastq -2 $(DATADIR)/$*.R2.fastq $(BT_OPTIONS_GLOBAL) -p $(THREADS) \
		| samtools sort -@ $(THREADS) -O sam - \
		| samtools view > $@

# Quality control single gzipped
$(DATADIR)/$(QCDIR)/%.sam: $(DATADIR)/%.U.fastq.gz | $(DATADIR)/$(QCDIR) Refs
	@echo -e $(red) "Running single read quality control for $*" $(end)
	bowtie2 -x refs/ctrl_capture -U $(DATADIR)/$*.U.fastq.gz $(BT_OPTIONS_GLOBAL) -p $(THREADS) \
		| samtools sort -@ $(THREADS) -O sam - \
		| samtools view > $@

# Quality control single
$(DATADIR)/$(QCDIR)/%.sam: $(DATADIR)/%.U.fastq | $(DATADIR)/$(QCDIR) Refs
	@echo -e $(red) "Running single read quality control for $*" $(end)
	bowtie2 -x refs/ctrl_capture -U $(DATADIR)/$*.U.fastq $(BT_OPTIONS_GLOBAL) -p $(THREADS) \
		| samtools sort -@ $(THREADS) -O sam - \
		| samtools view > $@

## Create dirs

# Create Visu dir
$(DATADIR)/$(VISUDIR):
	@echo -e $(red) "Creating $@" $(end)
	mkdir $@

# Create Blat dir
$(DATADIR)/$(BLATDIR):
	@echo -e $(red) "Creating $@" $(end)
	mkdir $@

# Create Fasta dir
$(DATADIR)/$(FASTADIR):
	@echo -e $(red) "Creating $@" $(end)
	mkdir $@

# Create Local alignment dir
$(DATADIR)/$(LOCALDIR):
	@echo -e $(red) "Creating $@" $(end)
	mkdir $@

# Create Global alignment dir
$(DATADIR)/$(GLOBALDIR):
	@echo -e $(red) "Creating $@" $(end)
	mkdir $@

# Create QC dir
$(DATADIR)/$(QCDIR):
	@echo -e $(red) "Creating $@" $(end)
	mkdir $@

## Bowtie2 references

# Build all the binary ref files
Refs: refs/HPVs.fasta $(patsubst refs/%.fasta, refs/%.1.bt2, $(wildcard refs/*.fasta))

# Compiling the references
refs/HPVs.fasta:
	@echo -e $(red) "Compiling the references" $(end)
	rm refs/HPVs.fasta -f
	cat refs/HPV*.fasta > refs/HPVs.fasta

# Build a binary ref file from fasta
refs/%.1.bt2 refs/%.2.bt2 refs/%.3.bt2 refs/%.4.bt2 refs/%.rev.1.bt2 refs/%.rev.2.bt2: refs/%.fasta
	@echo -e $(red) "Building binary reference for $*" $(end)
	bowtie2-build $< refs/$*

# Obtain hg19 ref from ucsc
hg19.2bit:
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

## Phony targets

.PHONY: All Refs QC_sam QC Global Genotyping Selection Local Profile Visu

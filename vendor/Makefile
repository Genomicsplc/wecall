SUBDIRS=samtools tabix

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) --directory=$@

clean:
	$(MAKE) --directory=samtools clean
	$(MAKE) --directory=tabix clean

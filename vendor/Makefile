SUBDIRS = samtools tabix
SUBDIRSCLEAN=$(addsuffix clean,$(SUBDIRS))

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean: $(SUBDIRSCLEAN)

%clean: %
	$(MAKE) -C $< clean

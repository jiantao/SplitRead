CURR_DIR := $(shell pwd)
SOURCES  := $(wildcard *.c)
DEPENDANTS = $(SOURCES:.c=.d)

dep: $(DEPENDANTS)

$(DEPENDANTS): $(SOURCES)

%.d: %.c
	@$(CC) -MM $(CFLAGS) $(INCLUDES) $< > $(OBJ_DIR)/$@; \
	echo '\t@$(CC) $(CFLAGS) $(INCLUDES) -c $(CURR_DIR)/$< -o $(OBJ_DIR)/$(*F).o\n' >> $(OBJ_DIR)/$@

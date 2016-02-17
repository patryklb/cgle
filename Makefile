CC=gcc
CXXFLAGS := -Wall -Wextra -pedantic-errors -std=c++0x
DEPS =
OBJ = main.o
LIBS=-lgsl -lgslcblas -lm -lstdc++ -lgomp

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
clean:
	rm -f $(OBJ)
	find . -depth -name "S*" -exec rm -rf '{}' \;

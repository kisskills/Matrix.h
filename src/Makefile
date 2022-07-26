CFLAGS=-Wall -Wextra -Werror -std=c11 #-fsanitize=address #-Wpedantic
TEST_FLAGS=-lcheck -lm -lpthread  #-lsubunit

all: s21_matrix.a gcov_report

rebuild: clean all

gcov_report: s21_matrix.c s21_helpers.c s21_test.c
	gcc $(CFLAGS) --coverage -o make_report.out s21_test.c s21_helpers.c s21_matrix.c $(TEST_FLAGS)
	./make_report.out
	lcov -t "make_report" -o make_report.info -c -d .
	genhtml -o report make_report.info

test: s21_matrix.o s21_helpers.o s21_test.o
	gcc -o test.out s21_matrix.o s21_helpers.o s21_test.o $(TEST_FLAGS)
	./test.out

s21_test.o: s21_test.c
	gcc $(CFLAGS) -c s21_test.c

s21_matrix.a: s21_matrix.o s21_helpers.o
	ar rcs s21_matrix.a s21_matrix.o s21_helpers.o
	ranlib s21_matrix.a

s21_matrix.o: s21_matrix.c
	gcc $(CFLAGS) -c s21_matrix.c

s21_helpers.o: s21_helpers.c
	gcc $(CFLAGS) -c s21_helpers.c

style: clean
	cppcheck *.h *.c
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extension=c *.c *.h

clean:
	rm -rf *.o *.a *.html *.gcda *.gcno *.css *.out *.info report

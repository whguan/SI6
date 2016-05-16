CC = icpc  -g
CFLAGS = #-O3 

C = $(CC) $(CFLAGS)
EXEC = ABwhguan

.PHONY : clean

$(EXEC):
	$(C) -c $(EXEC).cpp
	$(C) $(EXEC).o -o $(EXEC)
AB:
	./$(EXEC)  input.txt 1
clean:
	rm -f *.o
	rm -f $(EXEC)
	echo "clean obj"

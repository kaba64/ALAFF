CC := gcc 
LFLAGS = -lm
TARGET_EXEC := out
OBJS = main.o Gauss_Seidel.o Jacobi.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)

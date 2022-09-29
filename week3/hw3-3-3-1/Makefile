CC := gcc 
LFLAGS = -lflame -lopenblas -lm
TARGET_EXEC := out
OBJS = H_QR.o Housev.o axpy_inverse.o Norm_2_unb.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)

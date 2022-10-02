CC := gcc 
LFLAGS = -lflame -lopenblas -lm
TARGET_EXEC := out
OBJS = main.o LU_right_looking.o M_Rank1_Update.o G_VScale.o G_AXPY.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)

CC := gcc 
LFLAGS = -lflame -lopenblas -lm
TARGET_EXEC := out
OBJS = CGS_QR.o axpy_inverse.o dot.o Dot_Column_v_unb.o Mvmult_A_Tr_unb.o Mvmult_QR_unb.o Norm_2_unb.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)

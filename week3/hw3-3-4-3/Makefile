CC := gcc 
LFLAGS = -lflame -lopenblas -lm
TARGET_EXEC := out
OBJS = HQR.o Housev.o axpy_inverse.o Norm_2_unb.o Dot_Column_v_unb.o H_Rank1_Update.o Rank1-update.o Update_HQR_at12.o zero_vector.o Wt_HQR.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)

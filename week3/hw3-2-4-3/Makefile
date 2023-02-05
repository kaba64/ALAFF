CC := gcc 
LFLAGS = -lflame -lopenblas -lm
TARGET_EXEC := out
OBJS = MGS_QR.o axpy_inverse.o Dot_Column_v_unb.o Norm_2_unb.o Vmmult_q_Tr.o MGS_Rank1_Update.o Rank1_Update_Sub.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)

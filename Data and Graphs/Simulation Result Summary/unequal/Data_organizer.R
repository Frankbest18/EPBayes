# Data Table Path
EPB_Power_path = './Seperated Data/EPB/EPB_Power_summary.csv'
EPB_FDR_path = './Seperated Data/EPB/EPB_FDR_summary.csv'
Welch_Data_path = './Seperated Data/Welch Property Simulation Data/Welch Fails/Welch_Fails.csv'
Pooled_t_test_Data_path = './Seperated Data/Pooled_t_test.csv'
B_F_path = './Seperated Data/B_F.csv'

# Load Seperated Data Tables
EPB_Power_Data = read.csv(EPB_Power_path)
EPB_FDR_Data = read.csv(EPB_FDR_path)
Welch_Data = read.csv(Welch_Data_path)
Pooled_t_test_Data = read.csv(Pooled_t_test_Data_path)
B_F_Data = read.csv(B_F_path)

# Load Seperated Power Data
Power_1D_NPMLE = EPB_Power_Data[,c(1,2,5)]
Power_2D_NPMLE = EPB_Power_Data[, c(1,2,6)]
Power_Welch = Welch_Data[, c(2,3,4)]
Power_Pooled_t = Pooled_t_test_Data[, c(2,3,4)]
Power_B_F = B_F_Data[, c(1,2,3)]
Power_EV_NPMLE = EPB_Power_Data[, c(1,2,10)]
Power_EV2_NPMLE = EPB_Power_Data[, c(1,2,11)]

# Load Seperated FDR Data
FDR_1D_NPMLE = EPB_FDR_Data[,c(1,2,5)]
FDR_2D_NPMLE = EPB_FDR_Data[, c(1,2,6)]
FDR_Welch = Welch_Data[, c(2,3,5)]
FDR_Pooled_t = Pooled_t_test_Data[, c(2,3,5)]
FDR_B_F = B_F_Data[, c(1,2,4)]
FDR_EV_NPMLE = EPB_FDR_Data[, c(1,2,10)]
FDR_EV2_NPMLE = EPB_FDR_Data[, c(1,2,11)]

# Create Empty Data frames for United Data
Full_FDR_Data = data.frame('n1' = FDR_Welch$n1, 'n2' = FDR_Welch$n2, 'Welch' = FDR_Welch$fdr, 'Pooled_t' = FDR_Pooled_t$fdr, '1D_NPMLE' = rep(0, length(FDR_Welch$n1)), '2D_NPMLE' = rep(0, length(FDR_Welch$n1)), 'EV_NPMLE' = rep(0, length(FDR_Welch$n1)), 'EV2_NPMLE' = rep(0, length(FDR_Welch$n1)), 'B_F' = rep(0, length(FDR_Welch$n1)))
Full_Power_Data = data.frame('n1' = Power_Welch$n1, 'n2' = Power_Welch$n2, 'Welch' = Power_Welch$power, 'Pooled_t' = Power_Pooled_t$power, '1D_NPMLE' = rep(0, length(Power_Welch$n1)), '2D_NPMLE' = rep(0, length(Power_Welch$n1)), 'EV_NPMLE' = rep(0, length(Power_Welch$n1)), 'EV2_NPMLE' = rep(0, length(Power_Welch$n1)), 'B_F' = rep(0, length(FDR_Welch$n1)))

# Pair Separated Data based on (n1, n2) and Log them
for (i in c(1:nrow(Full_FDR_Data))) {
  n1 = Full_FDR_Data[i, 1]
  n2 = Full_FDR_Data[i, 2]
  
  for (j in c(1:nrow(EPB_FDR_Data))) {
    n1_EPB = EPB_FDR_Data[j, 1]
    n2_EPB = EPB_FDR_Data[j, 2]
    if (n1 == n1_EPB & n2 == n2_EPB) {
      Full_FDR_Data[i, 5] = EPB_FDR_Data[j, 5]
      Full_FDR_Data[i, 6] = EPB_FDR_Data[j, 6]
      Full_FDR_Data[i, 7] = EPB_FDR_Data[j, 10]
      Full_FDR_Data[i, 8] = EPB_FDR_Data[j, 11]
    }
  }
  
  for (j in c(1:nrow(B_F_Data))) {
    n1_B_F = B_F_Data[j, 1]
    n2_B_F = B_F_Data[j, 2]
    if (n1 == n1_B_F & n2 == n2_B_F) {
      Full_FDR_Data[i, 9] = B_F_Data[j, 4]
    }
  }
}

for (i in c(1:nrow(Full_Power_Data))) {
  n1 = Full_Power_Data[i, 1]
  n2 = Full_Power_Data[i, 2]
  
  for (j in c(1:nrow(EPB_Power_Data))) {
    n1_EPB = EPB_Power_Data[j, 1]
    n2_EPB = EPB_Power_Data[j, 2]
    if (n1 == n1_EPB & n2 == n2_EPB) {
      Full_Power_Data[i, 5] = EPB_Power_Data[j, 5]
      Full_Power_Data[i, 6] = EPB_Power_Data[j, 6]
      Full_Power_Data[i, 7] = EPB_Power_Data[j, 10]
      Full_Power_Data[i, 8] = EPB_Power_Data[j, 11]
    }
  }
  
  for (j in c(1:nrow(B_F_Data))) {
    n1_B_F = B_F_Data[j, 1]
    n2_B_F = B_F_Data[j, 2]
    if (n1 == n1_B_F & n2 == n2_B_F) {
      Full_Power_Data[i, 9] = B_F_Data[j, 3]
    }
  }
}

# Write United Simulation Data
write.csv(Full_Power_Data, './Full Data/Full_Power_Data_test.csv', row.names = FALSE)
write.csv(Full_FDR_Data, './Full Data/Full_FDR_Data_test.csv', row.names = FALSE)
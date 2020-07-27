import wannier_supercell as wsc

wsyst = wsc.wannier_supercell("linear_chain",[100,1,1])
wsyst.save_op_asCSR( op_list=["VX","SZ"] )

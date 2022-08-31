# MCPowerSEM: Structural Equation Modeling Power Analysis use Monte-Carlo Method #
#                                                                                #
# Author: Kazuha_Yae                                                             #
# Ref: https://journal.psych.ac.cn/xlkxjz/CN/10.3724/SP.J.1042.2022.02117        #
#                                                                                #
# Program Version: Alpha 0.0.1                                                   #
# R Version: ≥ 3.4                                                               #
# Depends: lavaan                                                                #
# License: GPL ≥ 2                                                               #
# Bug Report: sulfonamides@163.com                                               #
#================================================================================#
# 0. import packages
if(!require(lavaan))   install.packages('lavaan')
if(!require(mvtnorm))  install.packages('mvtnorm')
if(!require(progress)) install.packages('progress')
if(!require(scales))   install.packages('scales')

# 1. Set Class S3 MCPowerSEM
#--------------------------------------------------------------------------------
# 1.1 MCPowerSEM属性设定：
# 传入属性
#   generate.mod: 用以生成模拟数据的模型（以lavaan模型指定）
#   generate.cov: 用以生成模拟数据的协方差矩阵（直接以矩阵形式指定）
#   std.reg: 当以模型生成数据时，是否将输入模型的回归系数部分进行标准化处理
#   analysis.mod: 用以进行分析的模型
#   simulate.n: 样本容量
#   repeated: 重复次数
#   seed: 起始种子
# 模型内建属性
#   origin.mod: 原始模型
#   simulate.seed: 生成模拟数据用的种子
#   result.fitting: 模拟得到的模型拟合
#   result.paraest: 模拟得到的参数估计值
#   result.simcov: 模拟数据的协方差矩阵
#--------------------------------------------------------------------------------
# 1.2 MCPowerSEM函数
# 输入参数：
#   genMod: 用以生成虚拟数据的理论模型              
#   genCov: 用以生成虚拟数据的协方差矩阵            
#   anaMod: 实际进行分析时所执行的模型              
#   REP   : 模拟的重复次数                          
#   N     : 模拟样本量                              
#   std   : 是否需要根据“智能”还原标准化的协方差矩阵
#   seed  : 随机数种子                              
# 主要功能：New一个新的对象，并且绑定属性
# 输出结果：绑定并已完成基本处理的MCPowerSEM对象
# 主函数体：
MCPowerSEM = function(genMod = NULL, genCov = NULL, anaMod = NULL, REP = 1000L, N, std = FALSE, seed = NULL){
	# 简单检查一下输入问题
	# 没有穷尽可能（哪怕是比较可能）的输入端的错误，所以可能出现意外的报错
	# 特别地，没有检查输入的genMod是否合法，尤其是genMod与std之间的相合性（懒得盘逻辑了，有需要的可以自己加）
	if(is.null(genMod) & is.null(genCov)) stop('请指定用以生成模拟数据的基础模型')
	if(!is.null(genCov) & is.null(colnames(genCov)) & is.null(rownames(genCov))) stop('请指定变量名称')
	if((!is.null(genMod)) & (!is.null(genCov))) warning('您同时指定了两种生成数据的模型，这里默认使用了传入的协方差矩阵')
	if(is.null(anaMod)) stop('请指定用以分析数据的模型')
	if(!is.numeric(seed)) warning('指定的随机数种子并非数值型，没有成功指定随机数种子')
	if(length(seed) > 1) warning('输入的随机数种子长度大于1，这里使用了第一个元素指定随机数种子')
	
	# 生成用于产生数据的协方差矩阵
	if((is.null(genCov)) & (!std)){
		genCov = sem(genMod)@implied$cov[[1]]						# 调用模型内建的内隐协方差矩阵
		colnames(genCov) = sem(genMod)@Model@dimNames[[2]][[1]]		# 调用模型显变量残差矩阵的名字为协方差矩阵命名（这俩变量名一致）
	}else if((is.null(genCov)) & (std)){
		L = sem(genMod)@Model@GLIST$lambda
		T = sem(genMod)@Model@GLIST$theta
		B = sem(genMod)@Model@GLIST$beta
		P = sem(genMod)@Model@GLIST$psi
		if(is.null(B)){
			genCov = L %*% P %*% t(L) + T
			diag(genCov) = rep(1, ncol(genCov))						# 强行标准化显变量
		}else{
			I = diag(rep(1, ncol(B)))
			IB= solve(I - B)
			P = diag((rep(1, nrow(B)) %*% t(solve(I ** 2)))[1, ])	# 根据回归系数强行解出对应的潜变量协差阵
			genCov = L %*% IB %*% P %*% t(L %*% IB) + T
			diag(genCov) = rep(1, ncol(genCov))						# 强行标准化显变量
		}
		colnames(genCov) = sem(genMod)@Model@dimNames[[2]][[1]]
	}
	
	# 用generate covariate matrix跑一次结果
	fit0 = sem(model = anaMod, sample.cov = genCov, sample.nobs = N, likelihood = 'wishart')
	
	# YHY方法，参见Cheng和Wu（2017）的研究，该方法用来解决模型拟合参数的分布的期望值与非中心参数不等的情况
	cdimname = colnames(genCov)
	rdimname = rownames(genCov)
	Sigma= fit0@implied$cov[[1]]
	S    = genCov
	if(fitmeasures(fit0, 'chisq') != 0){	
		if(fitmeasures(fit0, 'chisq') <= fitmeasures(fit0, 'df')){
			genCov = Sigma
		}else{
			NCP = fitmeasures(fit0, 'chisq') - fitmeasures(fit0, 'df')
			YHY = function(a){
				weightedCov = a * S + (1-a) * Sigma
				fit = sem(model = anaMod, sample.cov = weightedCov, sample.nobs = N, likelihood = 'wishart')
				return(fitmeasures(fit, 'chisq') - NCP)
			}
			a = uniroot(YHY, c(0, 1))$root
			genCov = a * S + (1-a) * Sigma
		}
	}
	colnames(genCov) = cdimname
	rownames(genCov) = rdimname
	
	# 重新计算一次结果
	fit0 = sem(model = anaMod, sample.cov = genCov, sample.nobs = N, likelihood = 'wishart')
	
	# 初始化随机数生成器，产生每次循环中的随机数种子，随机数种子为随机的八位整数
	if(is.numeric(seed)) set.seed(seed[1])
	simseed = 1e8 + floor(runif(REP) * 9e8)
	
	# 实例化进度条
	progressbar = progress_bar$new(format = '   calculating [:bar] :percent   Will complete in :eta', total = REP, clear = FALSE)
	
	# 生成结果记录
	resfit = list()
	respar = list()
	rescov = list()
	
	# 循环采样并计算
	for(i in 1:REP){
		# 初始化随机数种子
		set.seed(simseed[i])
		# 生成随机样本
		Sample = rmvnorm(N, sigma = genCov)
		if(is.null(rownames(genCov))){colnames(Sample) = colnames(genCov)}else{colnames(Sample) = rownames(genCov)}
		# 建模，并返回结果
		fit = sem(anaMod, Sample, likelihood = 'wishart')
		resfit[[i]] = fitmeasures(fit)
		respar[[i]] = lavInspect(fit, what = 'list')
		rescov[[i]] = cov(Sample)
		# 进度条+1
		progressbar$tick()
	}
	
	# 绑定属性
	structure(list(
	generate.mod = genMod,
	generate.cov = genCov,
	std = std,
	analysis.mod = anaMod,
	simulate.n = N,
	repeated = REP,
	seed = seed[1], 
	simulate.seed = simseed, 
	result.fitting = resfit, 
	result.paraest = respar,
	result.simcov = rescov,
	origin.mod = fit0
	), class = 'MCPowerSEM')
}
#--------------------------------------------------------------------------------

# 2. 为MCPowerSEM类创建（伪）方法
# U2FsdGVkX19RDTlBTyjJjN/dSixgwVluel/9L9XCfM+DjWNhmcS164tJWOne8sfxlvazdJ8FftYuQNI6Df8sl6BTvZWvZXlvufWFB4U5mDEc4wIJa3UMiO9FB1aBER+j89DxP7X4oxQ=
#--------------------------------------------------------------------------------
# 2.1 Cover系列方法
# 此系列方法讨论的是参数覆盖/估计值分布的情况，包括三种：parCover从样本的角度分析样本协方差矩阵对总体协方差矩阵的参数覆盖情况；
# estCover从模型的角度分析基于模拟数据估计得到的模型参数对模型“真实”参数的覆盖情况；
# fitCover从模型的角度分析基于模拟数据得到的模型拟合值的分布
# 此系列方法的标准化输出：该系列方法的标准化输出为一矩阵，矩阵的每一行是一个特定的参数，每一列分别为“真实参数”（即用来生成样本的参数），
#   样本参数估计值的均值，样本参数估计值的标准差，样本参数估计值的中位数，以及样本参数估计值的各种百分位数
#--------------------------------------------------------------------------------
# 2.1.1 parCover
parCover = function(object, ...) UseMethod('parCover')
parCover.MCPowerSEM = function(object, Probs = c(0.005, 0.995, 0.025, 0.975)){
	genCov = object$generate.cov[lower.tri(object$generate.cov, TRUE)]
	result = matrix(0, nrow = object$repeated, ncol = length(genCov))
	for(i in 1:object$repeated) result[i, ] = object$result.simcov[[i]][lower.tri(object$result.simcov[[i]], TRUE)]
	output = matrix(0, nrow = length(genCov), ncol = 4 + length(Probs))
	output[, 1] = genCov
	output[, 2] = apply(result, 2, mean)
	output[, 3] = apply(result, 2, sd)
	output[, 4] = apply(result, 2, median)
	output[, 4 + 1:length(Probs)] = t(apply(result, 2, quantile, probs = Probs))
	colnames(output) = c('origin', 'mean', 'sd', 'median', percent(Probs))
	return(output)
}
# 2.1.2 estCover
estCover = function(object, ...) UseMethod('estCover')
estCover.MCPowerSEM = function(object, Probs = c(0.005, 0.995, 0.025, 0.975)){
	fit      = lavInspect(object$origin.mod, what = 'list')
	result.m = matrix(0, nrow = object$repeated, ncol = nrow(fit))
	result.sd= matrix(0, nrow = object$repeated, ncol = nrow(fit))
	for(i in 1:object$repeated){
		result.m[i, ] = object$result.paraest[[i]]$est
		result.sd[i,] = object$result.paraest[[i]]$se
	}
	output.m = matrix(0, nrow = nrow(fit), ncol = 3 + length(Probs))
	output.sd= output.m
	output.m[, 1] = apply(result.m, 2, mean)
	output.m[, 2] = apply(result.m, 2, sd)	
	output.m[, 3] = apply(result.m, 2, median)	
	output.sd[, 1] = apply(result.sd, 2, mean)
	output.sd[, 2] = apply(result.sd, 2, sd)	
	output.sd[, 3] = apply(result.sd, 2, median)
	output.m[, 3 + 1:length(Probs)] = t(apply(result.m, 2, quantile, probs = Probs))
	output.sd[, 3 + 1:length(Probs)] = t(apply(result.sd, 2, quantile, probs = Probs))
	colnames(output.m) = c('est.mean', 'est.sd', 'est.median', paste0('est.', percent(Probs)))
	colnames(output.sd) = c('se.mean', 'se.sd', 'se.median', paste0('se.', percent(Probs)))
	rownames(fit) = paste(fit[, 'lhs'], fit[, 'op'], fit[, 'rhs'])
	return(cbind(fit[, c('est', 'se')], output.m, output.sd))
}
# 2.1.3 fitCover
fitCover = function(object, ...) UseMethod('fitCover')
fitCover.MCPowerSEM = function(object, measures = 'all', Probs = c(0.005, 0.995, 0.025, 0.975)){
 	NAME = names(fitmeasures(object$origin.mod, measures))
 	result = matrix(0, nrow = object$repeated, ncol = length(NAME))
	if(measures == 'all'){
		for(i in 1:object$repeated) result[i, ] = object$result.fitting[[i]]
	}else{
		for(i in 1:object$repeated) result[i, ] = object$result.fitting[[i]][measures]
	}
 	output = matrix(0, nrow = length(NAME), ncol = 8)
	output[, 1] = fitmeasures(object$origin.mod, measures)
	output[, 2] = apply(result, 2, mean)
	output[, 3] = apply(result, 2, sd)
	output[, 4] = apply(result, 2, median)
	output[, 4 + 1:length(Probs)] = t(apply(result, 2, quantile, probs = Probs))
	colnames(output) = c('origin', 'mean', 'sd', 'median', percent(Probs))
	rownames(output) = NAME
	return(output)	
}
#--------------------------------------------------------------------------------
# 2.2 Test系列方法
# 此系列方法讨论的是基于某些决策标准时，REP次模拟中给出“拒绝”结果的几率，包括两种：estTest讨论的是对模型参数Wald检验的检验力分析（目前仅支持双侧检验的情形）
# fitTest讨论的是对模型拟合指数的检验力分析（需要研究者自行指定决策标准）
# 由于技术力原因，目前estTest仅支持Wald检验（即Delta法），fitTest需要研究者自行指定决策标准，通常情况下这些决策标准相对而言容易获得
#   但进行等效性检验时会麻烦一些，如果涉及等效性检验的情形，研究者可能需要先参考王阳等人（2020）的研究确定决策标准，而后进行分析
#--------------------------------------------------------------------------------
# 2.2.1 estTest
# 需要注意的是，这里暂时仅支持双侧检验
estTest = function(object, ...) UseMethod('estTest')
estTest.MCPowerSEM = function(object, alpha = c(0.05, 0.01)){
	fit      = lavInspect(object$origin.mod, what = 'list')
	result.m = matrix(0, nrow = object$repeated, ncol = nrow(fit))
	result.sd= matrix(0, nrow = object$repeated, ncol = nrow(fit))
	for(i in 1:object$repeated){
		result.m[i, ] = object$result.paraest[[i]]$est
		result.sd[i,] = object$result.paraest[[i]]$se
	}
	result.z = result.m / result.sd
	results  = list()
	for(i in 1:length(alpha)){
		results[[i]] = apply(result.z, c(1, 2), function(x) ifelse(abs(x) > qnorm(1-alpha[i]/2), 1, 0))
	}
	output   = fit[, c('est', 'se')]
	output   = cbind(output, z.value <- fit[, 'est'] / fit[, 'se'], p.value = 2*pnorm(-abs(z.value)))
	for(i in 1:length(alpha)) output = cbind(output, apply(results[[i]], 2, mean))
	colnames(output) = c('est', 'se', 'z.value', 'p.value', paste0('power.at.', alpha))
	rownames(output) = paste(fit[, 'lhs'], fit[, 'op'], fit[, 'rhs'])
	return(output)
}
# 2.2.2 fitTest
# 特别需要注意的是，对拟合指数的检验多数是单侧的，所以这里给出了左侧拒绝域和右侧拒绝域两个结果，需要研究者自行确定使用哪个值
fitTest = function(object, ...) UseMethod('fitTest')
fitTest.MCPowerSEM = function(object, measures = c('rmsea', 'rmsea', 'rmsea', 'rmsea.ci.upper', 'cfi', 'mfi'), crit = c(0.08, 0.1, 0.05, 0.08, 0.9, 0.85)){
	# 检测传入参数的合理性
	if(length(measures) != length(crit)) stop('measures与crit的长度需相等，请检查！')
 	NAME = names(fitmeasures(object$origin.mod, measures))
 	result = matrix(0, nrow = object$repeated, ncol = length(NAME))
	for(i in 1:object$repeated) result[i, ] = object$result.fitting[[i]][measures]
	powerL = numeric(length(NAME))
	powerR = numeric(length(NAME))
	for(i in 1:length(NAME)){
		powerL[i] = mean(result[, i] < crit[i])
		powerR[i] = mean(result[, i] > crit[i])
	}
	output = cbind(crit, powerL, powerR)
	rownames(output) = measures
	colnames(output) = c('crit', 'power.Left', 'power.Right')
	return(output)
}
#--------------------------------------------------------------------------------
# 3. 既有的summary方法
# 没啥好说的，就是打印个结果
#--------------------------------------------------------------------------------
summary.MCPowerSEM = function(object, Probs = c(0.005, 0.995, 0.025, 0.975), 
	cover.measures = 'all', test.measures = c('rmsea', 'rmsea', 'rmsea', 'rmsea.ci.upper', 'cfi', 'mfi'), 
	crit = c(0.08, 0.1, 0.05, 0.08, 0.9, 0.85), alpha = c(0.05, 0.01), digit = 5, detail = FALSE){
	cat('使用蒙特卡洛法对SEM模型的统计检验力进行分析\n')
	cat('模型参数：\n')
	cat('用以生成虚拟数据的协方差矩阵：\n')
	print(object$generate.cov)
	cat('用以进行分析的模型：\n')
	print(object$analysis.mod)
	cat('模拟样本量：', object$simulate.n, '\n')
	cat('重复次数：', object$repeated, '\n')
	cat('初始种子：', ifelse(is.numeric(object$seed), object$seed[1], '未设置'), '\n')
	cat('\n')
	cat('模型参数检验的统计检验力：\n')
	print(round(estTest(object, alpha), digit))
	cat('\n')
	cat('模型拟合检验的统计检验力：\n')
	print(round(fitTest(object, test.measures, crit), digit))
	if(detail){
		cat('\n')
		cat('参数覆盖情况\n')
		cat('协方差覆盖情况\n')
		print(round(parCover(object, Probs), digit))
		cat('\n')
		cat('模型参数估计值覆盖情况\n')
		print(round(estCover(object, Probs), digit))
		cat('\n')
		cat('模型拟合情况\n')
		print(round(fitCover(object, cover.measures, Probs), digit))
	}
}
#================================================================================

# 示例程序
# 模型来源于Sarah
# genMod = '
# Y  ~ 0.1 * C + - 0.2 * X + 0.2 * M1 + 0.2 * M2
# M1 ~ 0.1 * C + - 0.2 * X 
# M2 ~ 0.1 * C + - 0.2 * X 
# X  ~ 0.1 * C
# '
# 
# anaMod = '
# Y  ~ C + X + b1 * M1 + b2 * M2
# M1 ~ C + a1 * X 
# M2 ~ C + a2 * X 
# X  ~ C
# ab1 := a1 * b1
# ab2 := a2 * b2
# '
# 
# ana.test = MCPowerSEM(genMod = genMod, genCov = NULL, anaMod = anaMod, REP = 1000L, N = 401, std = TRUE, seed = 12345678)
# calculating [=========>---------------------------------------------------]   13%   Will complete in 41s
# 
# summary(ana.test)
# 
# 输出：
# 模式输出：
# 使用蒙特卡洛法对SEM模型的统计检验力进行分析
# 模型参数：
# 用以生成虚拟数据的协方差矩阵：
#             Y       M1       M2       X     C
# [1,]  1.00000  0.26496  0.26496 -0.2688 0.112
# [2,]  0.26496  1.00000  0.04640 -0.1920 0.080
# [3,]  0.26496  0.04640  1.00000 -0.1920 0.080
# [4,] -0.26880 -0.19200 -0.19200  1.0000 0.100
# [5,]  0.11200  0.08000  0.08000  0.1000 1.000
# 用以进行分析的模型：
# [1] "\nY  ~ C + X + b1 * M1 + b2 * M2\nM1 ~ C + a1 * X \nM2 ~ C + a2 * X \nX  ~ C\nab1 := a1 * b1\nab2 := a2 * b2\n"
# 模拟样本量： 401 
# 重复次数： 1000 
# 初始种子： 12345678 
# 
# 模型参数检验的统计检验力：
#                   est      se  z.value p.value power.at.0.05 power.at.0.01
# Y ~ C         0.09834 0.04611  2.13258 0.03296         0.577         0.338
# Y ~ X        -0.19826 0.04753 -4.17097 0.00003         0.989         0.948
# Y ~ M1        0.20932 0.04651  4.50044 0.00001         0.995         0.978
# Y ~ M2        0.20932 0.04651  4.50044 0.00001         0.990         0.961
# M1 ~ C        0.10020 0.04906  2.04236 0.04112         0.564         0.319
# M1 ~ X       -0.20202 0.04906 -4.11767 0.00004         0.987         0.936
# M2 ~ C        0.10020 0.04906  2.04236 0.04112         0.550         0.306
# M2 ~ X       -0.20202 0.04906 -4.11767 0.00004         0.977         0.927
# X ~ C         0.10000 0.04975  2.01008 0.04442         0.530         0.285
# Y ~~ Y        0.82477 0.05832 14.14214 0.00000         1.000         1.000
# M1 ~~ M1      0.95320 0.06740 14.14214 0.00000         1.000         1.000
# M2 ~~ M2      0.95320 0.06740 14.14214 0.00000         1.000         1.000
# X ~~ X        0.99000 0.07000 14.14214 0.00000         1.000         1.000
# C ~~ C        1.00000 0.00000      Inf 0.00000         1.000         1.000
# ab1 := a1*b1 -0.04229 0.01392 -3.03796 0.00238         0.957         0.744
# ab2 := a2*b2 -0.04229 0.01392 -3.03796 0.00238         0.942         0.726
# 
# 模型拟合检验的统计检验力：
#                crit power.Left power.Right
# rmsea          0.08      0.938       0.062
# rmsea          0.10      0.977       0.023
# rmsea          0.05      0.838       0.162
# rmsea.ci.upper 0.08      0.176       0.824
# cfi            0.90      0.000       1.000
# mfi            0.85      0.000       1.000

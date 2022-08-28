# MCPowerSEM: Structural Equation Modeling Power Analysis use Monte-Carlo Method #
#                                                                                #
# Author: Kazuha_Yae                                                             #
# Ref: https://journal.psych.ac.cn/xlkxjz/CN/10.3724/SP.J.1042.2022.02117        #
#                                                                                #
# Program Version: Alpha 0.0.1                                                   #
# R Version: ≥ 3.4                                                               #
# Depends: lavaan                                                                #
# License: GPL ≥ 2                                                               #
#================================================================================#
# import packages
if(!require(lavaan))   install.packages('lavaan')
if(!require(mvtnorm))  install.packages('mvtnorm')
if(!require(progress)) install.packages('progress')
if(!require(scales))   install.packages('scales')

# Set Class S3 MCPowerSEM
#--------------------------------------------------------------------------------#
# MCPowerSEM属性设定：
# 传入属性
# generate.mod: 用以生成模拟数据的模型（以lavaan模型指定）
# generate.cov: 用以生成模拟数据的协方差矩阵（直接以矩阵形式指定）
# std.reg: 当以模型生成数据时，是否将输入模型的回归系数部分进行标准化处理
# analysis.mod: 用以进行分析的模型
# simulate.n: 样本容量
# repeated: 重复次数
# seed: 起始种子
# 模型内建属性
# origin.mod: 原始模型
# simulate.seed: 生成模拟数据用的种子
# result.fitting: 模拟得到的模型拟合
# result.paraest: 模拟得到的参数估计值
# result.simcov: 模拟数据的协方差矩阵

# MCPowerSEM函数
# 输入参数：
# 主要功能：New一个新的对象，并且绑定属性
# 输出结果：绑定并已完成基本处理的MCPowerSEM对象
MCPowerSEM = function(genMod = NULL, genCov = NULL, anaMod = NULL, REP = 1000L, N, std = FALSE, seed = NULL){
	# 1.1 简单检查一下输入问题
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
		# PS. 这里如果是做参数检验的话已经可以用了，但要是做协差阵的话还需要使用YHY法“调节”一下，但这里还没写
	}else{
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
	
	# 初始化随机数生成器，产生每次循环中的随机数种子，随机数种子为随机的八位整数
	if(is.numeric(seed)) set.seed(seed[1])
	simseed = 1e8 + floor(runif(REP) * 9e8)
	
	# 实例化进度条
	progressbar = progress_bar$new(format = '   calculating [:bar] :percent   Will completion in :eta', total = REP, clear = FALSE)
	
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

# 创建（伪）方法
# U2FsdGVkX19RDTlBTyjJjN/dSixgwVluel/9L9XCfM+DjWNhmcS164tJWOne8sfxlvazdJ8FftYuQNI6Df8sl6BTvZWvZXlvufWFB4U5mDEc4wIJa3UMiO9FB1aBER+j89DxP7X4oxQ=
# A. 协方差矩阵覆盖情况
# 该方法给出模拟数据协差阵对生成数据协差阵的覆盖情况
parCover = function(object, ...) UseMethod('parCover')
parCover.MCPowerSEM = function(object, Probs = c(0.005, 0.995, 0.01, 0.99, 0.025, 0.975, 0.05, 0.95, 0.1, 0.9)){
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
# B. 参数估计情况
# 该方法给出对关键参数的估计，包括原始参数
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
	result.z = abs(result.m / result.sd)
	result.z95 = apply(result.z, c(1, 2), function(x) ifelse(x > qnorm(0.975), 1, 0))
	result.z99 = apply(result.z, c(1, 2), function(x) ifelse(x > qnorm(0.995), 1, 0))
	power95 = apply(result.z95, 2, mean)
	power99 = apply(result.z99, 2, mean)
	colnames(output.m) = c('est.mean', 'est.sd', 'est.median', paste0('est.', percent(Probs)))
	colnames(output.sd) = c('se.mean', 'se.sd', 'se.median', paste0('se.', percent(Probs)))
	return(cbind(fit[, c('lhs', 'op', 'rhs', 'est', 'se')], output.m, output.sd, power95, power99))
}
# C. 模型拟合情况
# 待完善
# fitCover = function(object, ...) UseMethod('fitCover')
# fitCover.MCPowerSEM = function(object, measures = 'all'){
# 	Rname = names(fitmeasures(object$origin.mod, measures)
# 	result = 
# 	
# }

# D. summary方法
# 待完善
summary.MCPowerSEM = function(object, Probs1 = c(0.005, 0.995, 0.01, 0.99, 0.025, 0.975, 0.05, 0.95, 0.1, 0.9), Probs2 = c(0.005, 0.995, 0.025, 0.975), measures = 'all'){
	# summary方法中尚不完善，目前只是一股脑把所有结果都打印出来，距离美化输出还有很长的路_(:3」∠)_
	print(parCover(object, Probs1))
	print(estCover(object, Probs2))
}

## 示例程序
#模型来源于Sarah
#genMod = '
#Y  ~ 0.1 * C + - 0.2 * X + 0.2 * M1 + 0.2 * M2
#M1 ~ 0.1 * C + - 0.2 * X 
#M2 ~ 0.1 * C + - 0.2 * X 
#X  ~ 0.1 * C
#'
#
#anaMod = '
#Y  ~ C + X + b1 * M1 + b2 * M2
#M1 ~ C + a1 * X 
#M2 ~ C + a2 * X 
#X  ~ C
#ab1 := a1 * b1
#ab2 := a2 * b2
#'
#
#ana.test = MCPowerSEM(genMod = genMod, genCov = NULL, anaMod = anaMod, REP = 1000L, N = 401, std = TRUE, seed = 12345678)
#
#模式输出：
#DBQ有点太丑了
#        origin        mean         sd      median        0.50%      99.50%         1.00%      99.00%         2.50%      97.50%         5.00%     95.00%      10.00%     90.00%
# [1,]  1.00000  1.00365510 0.07064933  1.00309051  0.841175508  1.19270372  0.8514962337  1.18076283  0.8678425067  1.14095641  0.8870470929  1.1196394  0.91203645  1.0915783
# [2,]  0.26496  0.26605273 0.05017445  0.26623358  0.137403439  0.39831410  0.1498117417  0.38612543  0.1738602689  0.36371278  0.1825759340  0.3461803  0.20117464  0.3286379
# [3,]  0.26496  0.26394568 0.05383408  0.26200908  0.130962660  0.42925691  0.1431458158  0.39527755  0.1640337753  0.37086896  0.1815897846  0.3501508  0.19789752  0.3315806
# [4,] -0.26880 -0.26884719 0.05166134 -0.26764151 -0.404411716 -0.14761749 -0.3939148784 -0.15058508 -0.3761891443 -0.16394471 -0.3555479005 -0.1886681 -0.33757131 -0.2073096
# [5,]  0.11200  0.11334688 0.05037082  0.11502773 -0.008902889  0.23074384 -0.0009831529  0.22698286  0.0119726161  0.20791563  0.0284057830  0.1958246  0.04595424  0.1765011
# [6,]  1.00000  1.00350638 0.07293035  1.00132159  0.809669148  1.19938624  0.8416097764  1.18152610  0.8676779238  1.15573992  0.8871555864  1.1239094  0.91171603  1.0973731
# [7,]  0.04640  0.04577332 0.04970678  0.04436772 -0.077018980  0.17378359 -0.0676025786  0.16213374 -0.0504430627  0.14896581 -0.0383322657  0.1295220 -0.01937146  0.1112688
# [8,] -0.19200 -0.19018001 0.05021733 -0.18990893 -0.318226251 -0.07094888 -0.3031258169 -0.07895927 -0.2845812862 -0.09416507 -0.2710841099 -0.1094379 -0.25540873 -0.1264587
# [9,]  0.08000  0.08348789 0.05049367  0.08441632 -0.046789812  0.20664502 -0.0360826660  0.20269775 -0.0145076855  0.17856264  0.0002650451  0.1655820  0.01802215  0.1484336
#[10,]  1.00000  0.99741155 0.07064397  0.99345040  0.826438887  1.20218373  0.8468810763  1.17987014  0.8653005321  1.14562136  0.8865645187  1.1207957  0.91068391  1.0868526
#[11,] -0.19200 -0.19155065 0.05335692 -0.19124461 -0.340985613 -0.05743802 -0.3205574167 -0.06916820 -0.3046301442 -0.08960675 -0.2801117712 -0.1038822 -0.25914557 -0.1235539
#[12,]  0.08000  0.08164227 0.05039787  0.08240314 -0.050608404  0.20227980 -0.0395799419  0.19666077 -0.0221301068  0.18235738 -0.0012254189  0.1633741  0.01656806  0.1462353
#[13,]  1.00000  0.99267033 0.07144908  0.98881378  0.824390487  1.20508323  0.8301179277  1.16892921  0.8582580401  1.13775622  0.8813433687  1.1126485  0.90294816  1.0845861
#[14,]  0.10000  0.10050565 0.05153953  0.09935756 -0.028555288  0.23706232 -0.0188111432  0.22572023 -0.0002433234  0.20464529  0.0223251991  0.1858816  0.03557371  0.1660760
#[15,]  1.00000  1.00249771 0.07151780  1.00131788  0.820151846  1.19074332  0.8393431007  1.16552427  0.8693876636  1.13980602  0.8875597704  1.1149838  0.91056614  1.0934211
#   lhs op   rhs         est         se    est.mean     est.sd  est.median     est.0.5%   est.99.5%      est.2.5%   est.97.5%    se.mean       se.sd  se.median     se.0.5%
#1    Y  ~     C  0.09833512 0.04611077  0.09868580 0.04574860  0.10098778 -0.009802546  0.20867787  0.0068337635  0.18603383 0.04617133 0.002342288 0.04610304 0.040468998
#2    Y  ~     X -0.19825626 0.04753244 -0.20035298 0.04702864 -0.19999447 -0.328209573 -0.08833271 -0.2935815299 -0.11032577 0.04779635 0.002366301 0.04767095 0.042234848
#3    Y  ~    M1  0.20931574 0.04651002  0.20970715 0.04537833  0.21225032  0.092058765  0.31694676  0.1197535473  0.29498598 0.04645953 0.002294544 0.04623810 0.041239697
#4    Y  ~    M2  0.20931574 0.04651002  0.20835276 0.04877340  0.20856015  0.075850250  0.33319009  0.1129307609  0.31047521 0.04661324 0.002231492 0.04661803 0.041318715
#5   M1  ~     C  0.10020202 0.04906180  0.10342442 0.04930387  0.10502718 -0.024680248  0.22136875  0.0098386568  0.19962616 0.04909718 0.002557189 0.04898140 0.043227944
#6   M1  ~     X -0.20202020 0.04906180 -0.20174806 0.04788109 -0.20211048 -0.327568131 -0.08362869 -0.2952924516 -0.10941553 0.04933657 0.002489941 0.04942175 0.043107000
#7   M2  ~     C  0.10020202 0.04906180  0.10169981 0.04895356  0.10055969 -0.020188186  0.21826813  0.0030590290  0.19741059 0.04892473 0.002405622 0.04880253 0.042843079
#8   M2  ~     X -0.20202020 0.04906180 -0.20306319 0.05150260 -0.20318432 -0.338826399 -0.07058143 -0.3061966406 -0.09941184 0.04917143 0.002506455 0.04911470 0.043416621
#9    X  ~     C  0.10000000 0.04974937  0.10004962 0.05069563  0.10107067 -0.032043506  0.22770972 -0.0002394369  0.20257887 0.04949676 0.002449203 0.04956511 0.043839684
#10   Y ~~     Y  0.82477459 0.05832037  0.81965905 0.05743941  0.82005375  0.690656441  0.98593081  0.7154189693  0.93694700 0.05795865 0.004061580 0.05798656 0.048836785
#11  M1 ~~    M1  0.95319596 0.06740113  0.95182244 0.06894067  0.94672094  0.784414609  1.13790635  0.8262229996  1.09276491 0.06730401 0.004874841 0.06694328 0.055466489
#12  M2 ~~    M2  0.95319596 0.06740113  0.94522274 0.06573843  0.94192240  0.785065312  1.13470340  0.8183053573  1.08825150 0.06683734 0.004648409 0.06660397 0.055512501
#13   X ~~     X  0.99000000 0.07000357  0.98003517 0.07035843  0.97765730  0.814233270  1.18963039  0.8463781411  1.11703052 0.06929895 0.004975092 0.06913081 0.057574987
#14   C ~~     C  1.00000000 0.00000000  1.00249771 0.07151780  1.00131788  0.820151846  1.19074332  0.8693876636  1.13980602 0.00000000 0.000000000 0.00000000 0.000000000
#15 ab1 := a1*b1 -0.04228601 0.01391923 -0.04225360 0.01352938 -0.04128599 -0.081839038 -0.01375365 -0.0697964887 -0.01857897 0.01415205 0.002222138 0.01418492 0.008624910
#16 ab2 := a2*b2 -0.04228601 0.01391923 -0.04237303 0.01522321 -0.04072567 -0.086293322 -0.01208903 -0.0789414096 -0.01697454 0.01416237 0.002442482 0.01409893 0.008130284
#     se.99.5%     se.2.5%   se.97.5% power95 power99
#1  0.05253891 0.041960506 0.05098144   0.577   0.338
#2  0.05425472 0.043373833 0.05262326   0.989   0.948
#3  0.05273971 0.042331165 0.05134280   0.995   0.978
#4  0.05242008 0.042367312 0.05103066   0.990   0.961
#5  0.05605161 0.044492733 0.05436957   0.564   0.319
#6  0.05575807 0.044336937 0.05456563   0.987   0.936
#7  0.05582106 0.044601731 0.05379084   0.550   0.306
#8  0.05655649 0.044514238 0.05429971   0.977   0.927
#9  0.05602675 0.045129066 0.05409101   0.530   0.285
#10 0.06971584 0.050587760 0.06625216   1.000   1.000
#11 0.08046213 0.058422789 0.07727015   1.000   1.000
#12 0.08023565 0.057862927 0.07695100   1.000   1.000
#13 0.08411957 0.059847972 0.07898599   1.000   1.000
#14 0.00000000 0.000000000 0.00000000   1.000   1.000
#15 0.01958879 0.009812065 0.01843948   0.957   0.744
#16 0.02062868 0.009525844 0.01935489   0.942   0.726

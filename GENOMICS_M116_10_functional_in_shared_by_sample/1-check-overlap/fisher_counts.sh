# for voi differences in m116 and 1000g

pop=filter(ref_overtable, CATEGORY==category)$
m=filter(overtable, CATEGORY==ref_category)$IN_CAT
k=filter(overtable, CATEGORY==category)$TOT_CAT
q=filter(overtable, CATEGORY==category)$IN_CAT
cont_mat <- matrix(c(q, m-q, k-q, n-(k-q)), 2, 2)

#f_test <- fisher.test(cont_mat, alternative=alternative, conf.int = TRUE, conf.level = 0.95)
f_test <- exact2x2::exact2x2(cont_mat)

# voi in both joint
pop=






## minimal object accepted by clis: .bic_design reads <par>.x directly
mkfit <- function(fam, y, Xmu, Xsi, Xnu, b, dl, g) {
  structure(list(family = c(fam, "x"), y = y,
    mu.fv = as.vector(plogis(Xmu %*% b)), sigma.fv = as.vector(exp(Xsi %*% dl)),
    nu.fv = as.vector(plogis(Xnu %*% g)),
    mu.link="logit", sigma.link="log", nu.link="logit",
    mu.x=Xmu, sigma.x=Xsi, nu.x=Xnu,
    mu.coefficients=b, sigma.coefficients=dl, nu.coefficients=g),
    class = c("gamlss","gam","glm","lm"))
}

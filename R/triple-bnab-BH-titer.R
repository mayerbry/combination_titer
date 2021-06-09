# Bryan Mayer
# 5-20-2021

# BH titer: common factor applied to each individual Ab concentration that reduces
# neutralization to 50%


#' Theoretical calculation for experimental (pooled sera) neutralization titer
#' under BH interaction
#'
#' @param ID50_1 ratio of concentration to IC50 for first Ab
#' @param ID50_2 ratio of concentration to IC50 for second Ab
#' @param ID50_3 ratio of concentration to IC50 for third Ab
#' @param titer_target goal combination titer endpoint: 0.5 = ID50, 0.8 = ID80
#'
#' This function is vectorized for use over PK
#'
#' @return
#' @export
#'
#' @examples
calc_3bnab_BHtiter = function(ID50_1, ID50_2, ID50_3, titer_target = 0.5, NaN_as_NA = F){
  stopifnot(titer_target > 0 & titer_target < 1)
  # if a row is all zeros, this is tracked then 0 is returned
  all_zeros = which(ID50_1 < sqrt(.Machine$double.eps) &
                  ID50_2 < sqrt(.Machine$double.eps) &
                  ID50_3 < sqrt(.Machine$double.eps))
  # any zeros has a closed-form solution
  any_zeros = which(ID50_1 < sqrt(.Machine$double.eps) |
                      ID50_2 < sqrt(.Machine$double.eps) |
                      ID50_3 < sqrt(.Machine$double.eps))
  # managing numerical precision
  ID50_1 = pmax(ID50_1, sqrt(.Machine$double.eps))
  ID50_2 = pmax(ID50_2, sqrt(.Machine$double.eps))
  ID50_3 = pmax(ID50_3, sqrt(.Machine$double.eps))

  id50_dat = data.frame(a = ID50_1, b = ID50_2, c = ID50_3)
  id50_dat$max_row = pmax(id50_dat$a, id50_dat$b, id50_dat$c)
  id50_dat$min_row = pmin(id50_dat$a, id50_dat$b, id50_dat$c)

  id50_dat$ID = NA_real_

  # if there are any zeros, the information only exists in at most two of the input titers
  if(length(any_zeros) > 0){

    #browser()
    id50_dat$ID[any_zeros] = .calc_2bnab_BHIDxx(ID50a = id50_dat$min_row[any_zeros],
                                                ID50b = id50_dat$max_row[any_zeros],
                                                titer_target = titer_target,
                                                NaN_as_NA = NaN_as_NA)
    id50_dat$ID[-any_zeros] = .calc_3bnab_BHIDxx(id50_dat[-any_zeros,],
                                                 titer_target = titer_target,
                                                 NaN_as_NA = NaN_as_NA)
  } else{
    id50_dat$ID = .calc_3bnab_BHIDxx(id50_dat,
      titer_target = titer_target,
      NaN_as_NA = NaN_as_NA)
  }

  id50_dat$ID[all_zeros] = 0
  id50_dat$ID
}


calc_2bnab_BHtiter = function(ID50_1, ID50_2, titer_target = 0.5, NaN_as_NA = F){
  zeros = which(ID50_1 < sqrt(.Machine$double.eps) &
                  ID50_2 < sqrt(.Machine$double.eps))  # managing numerical precision
  ID50_1 = pmax(ID50_1, sqrt(.Machine$double.eps))
  ID50_2 = pmax(ID50_2, sqrt(.Machine$double.eps))
  out = .calc_2bnab_BHIDxx(ID50_1, ID50_2, titer_target, NaN_as_NA)
  if(length(zeros > 0)) out[zeros] = 0
  out
}

.calc_2bnab_BHIDxx = function(ID50a, ID50b, titer_target = 0.5, NaN_as_NA = F){

  x = (-(ID50a + ID50b) + sqrt((ID50a + ID50b)^2 + 4*ID50b*ID50a * (titer_target/(1 - titer_target))))/(2*ID50a*ID50b)
  if(any(is.nan(x))) {
    if(!NaN_as_NA) stop(paste("BH ID50 nAn, input:", ID50a, ID50b))
    x <- NA_real_
  }
  out = 1/x

  # this handles some cluster issues with high values
  infinites = which(is.infinite(out))
  if(length(infinites > 0)) out[infinites] = pmax(ID50a, ID50b)[infinites]

  out

}

.BH_cubic_fun = function(x, ID50a, ID50b, ID50c, titer_target = 0.5){
  termA = ID50a * ID50b * ID50c
  termB =  ID50a * ID50b +  ID50a * ID50c +  ID50b * ID50c
  termC =  ID50a + ID50b + ID50c
  termD = titer_target/(1 - titer_target)
  x^3 * termA + x^2 * termB + x * termC - termD
}

.calc_3bnab_BHIDxx = function(ID50dat, titer_target = 0.5, NaN_as_NA = F){

  x = purrr::map_dbl(1:nrow(ID50dat), function(i){

    res = uniroot(.BH_cubic_fun, interval = c(0, 1),
                  ID50a = ID50dat$a[i], ID50b = ID50dat$b[i], ID50c = ID50dat$c[i],
                  titer_target = titer_target, tol = .Machine$double.eps ^ 0.5)
    res$root

  })


  if(any(is.nan(x))) {
    if(!NaN_as_NA) stop(paste("BH ID50 nAn, input:", ID50a, ID50b))
    x <- NA_real_
  }
  out = 1/x

  # this handles some cluster issues with high values
  infinites = which(is.infinite(out))
  if(length(infinites > 0)) out[infinites] = pmax(ID50a, ID50b, ID50c)[infinites]

  out

}

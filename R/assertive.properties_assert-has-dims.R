#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_has_cols <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                       
  assert_engine(
    has_cols, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_has_dims <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                                
  assert_engine(
    has_dims, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_has_rows <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                           
  assert_engine(
    has_rows, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

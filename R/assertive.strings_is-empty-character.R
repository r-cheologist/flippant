#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base call_and_name
#' @importFrom assertive.base set_cause
#' @importFrom assertive.base get_name_in_parent
is_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- !nzchar(x)
      set_cause(ok, ifelse(is.na(x), "missing", "nonempty"))
    },
    x
  )
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base call_and_name
#' @importFrom assertive.base set_cause
#' @importFrom assertive.base get_name_in_parent
is_non_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  call_and_name(
    function(x)
    {
      ok <- nzchar(x)
      set_cause(ok, "empty")
    },
    x
  )
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base is_na
#' @importFrom assertive.base set_cause
#' @importFrom assertive.base get_name_in_parent
is_missing_or_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  ok <- !nzchar(x) | is_na(x)
  set_cause(ok, "nonempty")
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base is_na
#' @importFrom assertive.base set_cause
#' @importFrom assertive.base get_name_in_parent
is_non_missing_nor_empty_character <- function(x, .xname = get_name_in_parent(x))
{ 
  x <- coerce_to(x, "character", .xname)
  ok <- nzchar(x) & !is_na(x)
  set_cause(ok, ifelse(is.na(x), "missing", "empty"))
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
is_not_missing_nor_empty_character <- function(x)
{
  .Deprecated("is_non_missing_nor_empty_character")
  is_non_missing_nor_empty_character(x)
}



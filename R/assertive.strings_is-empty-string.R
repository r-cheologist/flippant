#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base false
is_an_empty_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_a_string(x))) return(ok)
  if(!is_empty_character(x)) 
  {
    return(
      false(
        gettext("%s contains characters."), 
        .xname
      )
    )
  }
  TRUE
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base false
is_a_non_empty_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_a_string(x))) return(ok)
  if(!is_non_empty_character(x))
  {
    return(
      false(
        gettext("%s has no characters."), 
        .xname
      )
    )
  }
  TRUE
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base false
is_a_missing_or_empty_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_a_string(x, .xname))) return(ok)
  if(!is_missing_or_empty_character(x)) 
  {
    return(
      false(
        gettext("%s is not missing and contains characters."),
        .xname
      )
    )
  }
  TRUE
}

#' @author Richard Cotton <richierocks@gmail.com>
#' @noRd
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base false
is_a_non_missing_nor_empty_string <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_a_string(x))) return(ok)
  if(!is_non_missing_nor_empty_character(x)) 
  {
    return(
      false(
        gettext("%s is missing or has no characters."),
        .xname
      )
    )
  }
  TRUE
}

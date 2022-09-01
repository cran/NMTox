#' NM-400 in vitro dataset
#'
#' This dataset contains the result of genetic toxicity in vitro study of NM-400
#' (Multi-walled carbon nanotubes) with associated controls and variables
#' related to the experiments.
#'
#' @usage data(nm400)
#'
#' @format A data frame with columns:
#' \itemize{
#'  \item{name: }{Project-assigned name of the nanomaterial}
#'  \item{publicname: }{A widely accepted unique identifier}
#'  \item{supplier: }{Supplier/project where the data is originated from}
#'  \item{experimentID: }{Identifier of the experiment}
#'  \item{method: }{Method/assay used in the experiment}
#'  \item{studyprovider: }{Study provider}
#'  \item{endpoint: }{Toxicity endpoint measure}
#'  \item{value: }{Endpoint value}
#'  \item{unit: }{Unit of the endpoint}
#'  \item{celltype: }{Type of the cell used in the experiment}
#'  \item{treatment: }{Indicator of the treatment}
#'  \item{exptimeunit: }{Unit of measurement of the exposure time}
#'  \item{exptime: }{Exposure time}
#'  \item{concentration_unit: }{Unit of measurement of the concentration in
#'  variable \emph{concentration}}
#'  \item{concentration: }{Concentration of the nanomaterial}
#'  \item{concentration_ml_unit: }{Unit of measurement of the concentration
#'   in variable \emph{concentration_ml}}
#'  \item{concentration_ml: }{Concentration of the nanomaterial in amount per
#'  ml}
#' }
#'
#'
#' @source This dataset was obtained from
#'   \url{https://www.anses.fr/en/content/nanogenotox-project} (NanoGenotox
#'   project) and it was extracted from eNanoMapper database
#'   \url{https://search.data.enanomapper.net/}
#'
#'   The NANOGENOTOX Joint Action received funding from the European Union, in
#'   the framework of the Health Programme under Grant Agreement n2009 21.
#'
#'   Supported by European Union's Horizon 2020 research and innovation
#'   programme under grant agreement No 814426 - NanoInformaTIX
#'   \url{https://www.nanoinformatix.eu/}
#'
#' @examples
#' data(nm400)
#'
#'
#'
"nm400"

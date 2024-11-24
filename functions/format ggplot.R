# add title, subtitle and caption to a ggplot
format.ggplot <- function(plot.name,
                          title = experiment.title,
                          subtitle = experiment.subtitle,
                          caption = "") {
  formatted.plot <- plot.name +
    labs(
      title = title,
      subtitle = subtitle,
      caption = paste(caption, "Generated:", date())
    )
  return(formatted.plot)
}
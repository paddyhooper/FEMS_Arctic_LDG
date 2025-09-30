#Figure 6. Functional relative abundance plot
#Author: Patrick M. Hooper

func_RA_table <- read.table("~/func_ra_table.txt", header = TRUE)
func_RA_table

#Now we need to order by Latitude
region_order <- c("Kuujjuarapik", "Umiujaq", "Cambridge_Bay", "Bylot_Island", "Resolute", "Ellesmere_Island", "Ice_Shelves")

ggplot(func_RA_table, aes(fill=Nutrition, x = Region, y = Value)) + 
  geom_bar(position="stack", stat="identity")+ scale_x_discrete(limits=region_order)+ scale_fill_manual(values = function_palette)+ font_size + theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1))+labs(y= "Relative Abundance (%)", x = "Sampling Region")

#use ggsave to save to your wd as a pdf and jpeg
ggsave("function_RA.jpg", width = 297, height = 210, units = c("mm"))
ggsave("function_RA.pdf", width = 297, height = 210, units = c("mm"))

#### END OF SCRIPT ####

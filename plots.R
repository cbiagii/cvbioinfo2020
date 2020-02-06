library(KEGGREST)
listDatabases()


ko <- keggList("ko")
ko <- data.frame(ko = names(ko), 
                 description = ko)

pathway <- keggList("pathway")
pathway <- data.frame(pathway = names(pathway), 
                      description = pathway)


#Pathway
sig_path <- subset(dataset_path_scale, abs(mean_all) > 1)
sig_path$Description <- pathway$description[unique(grep(paste(gsub("ko", "", sig_path$Pathway), collapse = "|"), pathway$pathway))]

ggplot(sig_path, aes(x=reorder(Description, mean_all), y=mean_all, label=mean_all)) + 
  geom_bar(stat='identity', aes(fill=type), width=.5)  +
  scale_fill_manual(name="Count", 
                    labels = c("Above Average", "Below Average"), 
                    values = c("above"="#00ba38", "below"="#f8766d")) + 
  coord_flip() + labs(x="", y="Mean Counts", title="Pathways")


#Orthology
sig_path <- subset(dataset_path_scale, abs(mean_all) > 1)
sig_path$Description <- ko$description[unique(grep(paste(gsub("ko", "", sig_path$Pathway), collapse = "|"), ko$ko))]
ggplot(sig_path, aes(x=reorder(Description, mean_all), y=mean_all, label=mean_all)) + 
  geom_bar(stat='identity', aes(fill=type), width=.5)  +
  scale_fill_manual(name="Count", 
                    labels = c("Above Average", "Below Average"), 
                    values = c("above"="#00ba38", "below"="#f8766d")) + 
  coord_flip() + labs(x="", y="Mean Counts", title="Orthology")

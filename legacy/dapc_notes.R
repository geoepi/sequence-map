
merge polygons with PCs.  Making a separate copy of the polygons.
```{r}
vn_provinces_pc <- vn_provinces %>%
  left_join(unsup_pcs, 
            by = c("prov_eng" = "location"))
```


### View PC1  
Any spatial patterns in how the values are distributed? Gray areas are locations without samples. 
```{r}
ggplot() +
  geom_sf(data = vn_provinces_pc, aes(fill = PC_1), color = "black", size = 0.5) +
  scale_fill_viridis_c(name = "PC-1", 
                       option = "turbo",
                       na.value = "gray90") +
  coord_sf() +
  labs(title = " ",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
    legend.position = "right",
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(2, "line"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x =  element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 22, face = "bold"))
```

### View all PCs 
Plotting together may not be optimal because some PCs may have much higher or lower value ranges than others, causing differences to be washed out or not viable.
```{r}
vn_provinces_long <- vn_provinces_pc %>%
  pivot_longer(
    cols = starts_with("PC_"),         
    names_to = "Principal_Component", 
    values_to = "PC_Value"    
  )

ggplot(data = vn_provinces_long) +
  geom_sf(aes(fill = PC_Value), color = "black", size = 0.5) +
  scale_fill_viridis_c(name = "PC Value", 
                       option = "turbo",
                       na.value = "gray90") +
  coord_sf() +
  labs(title = "Principal Components across Provinces",
       x = "Longitude",
       y = "Latitude") +
  facet_wrap(~Principal_Component, ncol = 4) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
    legend.position = "right",
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(2, "line"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x =  element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 22, face = "bold"),
    strip.text = element_text(size = 16, face = "bold") 
  )
```

## DF Maps

Add province to PCs
```{r}
unsup_dfs <- left_join(unsup_dfs, seq_meta, by = "accession")

head(unsup_dfs)
```

merge polygons with PCs.  Making a separate copy of the polygons.
```{r}
vn_provinces_df <- vn_provinces %>%
  left_join(unsup_dfs, 
            by = c("prov_eng" = "location"))
```


### View LD1  
Any spatial patterns in how the values are distributed? Gray areas are locations without samples. 
```{r}
ggplot() +
  geom_sf(data = vn_provinces_df, aes(fill = LD1), color = "black", size = 0.5) +
  scale_fill_viridis_c(name = "LD-1", 
                       option = "turbo",
                       na.value = "gray90") +
  coord_sf() +
  labs(title = " ",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
    legend.position = "right",
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(2, "line"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x =  element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 22, face = "bold"))
```

### View all PCs 
Plotting together may not be optimal because some PCs may have much higher or lower value ranges than others, causing differences to be washed out or not viable.
```{r}
vn_provinces_long <- vn_provinces_pc %>%
  pivot_longer(
    cols = starts_with("PC_"),         
    names_to = "Principal_Component", 
    values_to = "PC_Value"    
  )

ggplot(data = vn_provinces_long) +
  geom_sf(aes(fill = PC_Value), color = "black", size = 0.5) +
  scale_fill_viridis_c(name = "PC Value", 
                       option = "turbo",
                       na.value = "gray90") +
  coord_sf() +
  labs(title = "Principal Components across Provinces",
       x = "Longitude",
       y = "Latitude") +
  facet_wrap(~Principal_Component, ncol = 4) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"),
    legend.position = "right",
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(2, "line"),
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 22, face = "bold"),
    axis.text.x =  element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 22, face = "bold"),
    strip.text = element_text(size = 16, face = "bold") 
  )
```




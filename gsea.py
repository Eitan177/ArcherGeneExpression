import os
import pandas as pd
import numpy as np
import re
import gseapy as gp
import streamlit as st
import seaborn as sns
import matplotlib.pyplot as plt



with st.form(key='file upload'):
    # Add form elements
    uploaded_file = st.file_uploader("Upload archer file")
    selected_option = st.selectbox('Select an option', ['minimum', 'maximum', 'mean'])
    samplenorm = st.checkbox('normalize within sample')

    # Add a form submit button
    submitted = st.form_submit_button('Submit')

# Process form submission
if submitted:
    if uploaded_file is not None:

        
        expressions=pd.read_table(uploaded_file,header=0,sep='\t')
        colnamesuse=expressions.columns
        expressions.columns=expressions.iloc[0]
        expressions.drop(0,axis=0,inplace=True)
        
        
        expressions['gene_abbrev']=[re.sub('_.+','',x) for x in expressions['GSP2']]
        expressions=expressions.iloc[:,np.where([str(x)!='nan' for x in expressions.columns])[0]] 
       
        e=expressions.iloc[:,np.arange(1,expressions.shape[1]-1,6)]

        #colnamesuse=pd.read_table(uploaded_file,sep='\t').columns  
        colnamesuse=colnamesuse[np.arange(1,len(colnamesuse)-1,6)]
        #breakpoint()
        e.columns=colnamesuse[0:e.shape[1]]
        e['gene_abbrev']=expressions['gene_abbrev']

       
        if selected_option=='maximum':
            gseaready=e.groupby('gene_abbrev')[colnamesuse].apply(lambda x: x.apply(lambda x: max([float(zz) for zz in x])))
        elif selected_option=='minimum':
            
            gseaready=e.groupby('gene_abbrev')[colnamesuse].apply(lambda x: x.apply(lambda x: min([float(zz) for zz in x])))

        elif selected_option=='mean':
            gseaready=e.groupby('gene_abbrev')[colnamesuse].apply(lambda x: x.apply(lambda x: np.mean([float(zz) for zz in x])) )
            
        gseaready.columns=[re.sub('-.+','',x) for x in gseaready.columns]
        try:
            gseaready.drop('gene_abbrev',axis=1,inplace=True)
        except:
            st.write('no gene abbrev column name')    
        if samplenorm:
            gseaready=gseaready.apply(lambda x: np.array([float(zz) for zz in x])/sum([float(zz) for zz in x]),axis=0)

        heatmap = sns.heatmap(gseaready, annot=False, cmap="YlGnBu")
        ftabs=st.tabs(['Gene Level no reordering','Gene Level with reordering','Pathway no reordering','Pathway with reordering'])
        # Display the heatmap in the Streamlit app
        with ftabs[0]:
            plt.figure(figsize=(10, 7))
            st.title('Gene Level Heatmap without reordering')
            st.pyplot(heatmap.get_figure())
        with ftabs[1]: 
            plt.figure(figsize=(10, 7))
            clustermap=sns.clustermap(gseaready, annot=False, cmap="YlGnBu")
            st.title('Gene Level Heatmap with reordering')
            st.pyplot(clustermap)

# Show the graph

        #st.plotly_chart(fig, theme="streamlit",height=1000,use_container_width=True)
        pthwylist=[]
        #tas=st.tabs(gseaready.columns.tolist())
        m=0
        
        for cc in gseaready.columns:

            gg=pd.DataFrame(gseaready[cc])
            tt=gp.ssgsea(data=gg,gene_sets="c2.all.v2023.2.Hs.symbols.gmt")
            #with tas[m]:
            #    st.write(tt.res2d.sort_values('Name'))
            pthwylist.append(tt.res2d.sort_values('Term'))
                
            m+=1
        pthwy=pd.concat(pthwylist,axis=1)   
        enriched=pd.DataFrame(pthwy[['NES']])
        enriched.columns=pthwy[['Name']].iloc[0].to_list()
    
        enriched.index=pthwy[['Term']].iloc[:,0].apply(lambda x:x[0:50])
        enriched=enriched.apply(lambda x: [float(zz) for zz in x])


        # Display the heatmap in the Streamlit app
        with ftabs[2]:
            plt.figure(figsize=(10, 7))
            heatmap = sns.heatmap(enriched, annot=False, cmap="YlGnBu")
            st.title('Pathway Heatmap without reordering')
            st.pyplot(heatmap.get_figure())
        with ftabs[3]:    
            plt.figure(figsize=(10, 7))
            clustermap=sns.clustermap(enriched, annot=False, cmap="YlGnBu")
            st.title('Pathway Heatmap with reordering')
            st.pyplot(clustermap)


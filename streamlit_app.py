from model import *
import streamlit as st

tab1, tab2, tab3 = st.tabs(["model", "example", "how it works"])

with tab1:
  st.title("3D BU Statistical Model")
  st.image("3D Blender.png")
  st.write("Hi and welcome, this tool creates a model of a blender depending on your inputs. It's recommended to look at the :blue-background[example] tab if this is your first time.")
  st.write("If you have any feedback, please reach out to either **Isaiah Gonzales** or **Rajarshi Sengupta**.")
           
  st.sidebar.write("**Input simulation parameters here**")
  model_type = st.sidebar.selectbox("Would you like to simulate one blender or multiple blenders", ["-","single run","multiple runs"], help ="**Single:** This will simulate one blender, and return *individual* assays of BU samples extracted. **Multiple:** This will simulate multiple blenders, and return *mean* assays, this can be useful to understand probabilities. ") 
  
  if model_type != "-":
    distribution = st.sidebar.selectbox("Please choose how you wish the powder to be distributed in the blender", ["unmixed", "random", "uniform", "poor"])
    blenderSize = st.sidebar.number_input("Size of blender (mL)", min_value=5, max_value=100000)
    advanced_options = st.sidebar.expander("Advanced Options")
    with advanced_options:
      thiefSize = st.slider("Size of sample thief (mL)", min_value = 1, max_value = 10, step =1)
      percentPurityOfDS = st.slider("Purity of DS (%)", min_value= 0, max_value = 110, value = 100, step=10)
      DL = st.slider("Blend drug load (%)", min_value = 0, max_value = 100, value = 20, step=10)
      fillRatio = st.slider("Fill Ratio", min_value=0.4, max_value=0.9, step= 0.1, value=0.5, help="What ratio of the blenders total volume is filled with powder?")
      particleSize = st.slider("Particle Size", min_value = 50, max_value=200, step = 10, value = 100)
  
    if distribution == "poor":
      with advanced_options:
        percentClumps = st.slider("What percent of DS particles would you like clumped?", min_value = 1,value=50, step =10)
        sizeClumps = st.number_input("Size of DS clumps (microns)", min_value = 100,value=1000, max_value =10000)
    
  if model_type == "multiple runs":
    numLoops = int(st.sidebar.number_input("How many simulations would you like to perform and average", min_value = 1, max_value = 500))

  runSimButton = st.sidebar.button("Get my sampling results.")
  if st.session_state.get('button') != True:
    st.session_state['button'] = runSimButton
  
  if st.runSimButton:
    if model_type == "single run":
      if distribution == "poor":
        blender3D(
                  thiefSize=thiefSize, 
                  percentPurityOfDS=percentPurityOfDS,
                  DL=DL, 
                  blenderSize=blenderSize, 
                  fillRatio=fillRatio/100,
                  distribution= distribution,
                  clumpiness=percentClumps/10,
                  clumpSize=sizeClumps,
                  particleSize = particleSize)
      else:
        blender3D(
                  thiefSize=thiefSize, 
                  percentPurityOfDS=percentPurityOfDS,
                  DL=DL, 
                  blenderSize=blenderSize, 
                  fillRatio=fillRatio,
                  distribution= distribution,
                  particleSize = particleSize)

  if st.session_state['button'] == True:
    if st.button("Visualize my model (this can take a while)"):
      displayBlender(placeholderaxes, blender,filledspace, top, middle, bottom, TopSamplingArray, MidSamplingArray, BotSamplingArray, particleSize, distribution, percentPurityOfDS)
  
    # if model_type == "multiple runs":
    #   results = []
    #   meanResults = []
    #   i = 0
    #   progbar = st.progress(0.0)
    #   if distribution == "poor":
    #     while i < numLoops:
    #       progbar.progress(i/numLoops, text = "simulation " + str(i) + " of " + str(numLoops))
    #       results.append(blender3D(thiefSize=thiefSize, 
    #                 percentPurityOfDS=percentPurityOfDS,
    #                 DL=DL, 
    #                 blenderSize=blenderSize, 
    #                 fillRatio=(fillRatio/100),
    #                 distribution= distribution,
    #                 clumpiness=percentClumps/10,
    #                 clumpSize=sizeClumps)
    #       i += 1
    #     progbar.empty()
    #   else:
    #     while i < numLoops:
    #       progbar.progress(i/numLoops, text = "simulation " + str(i) + " of " + str(numLoops))
    #       results.append(blender3D(thiefSize=thiefSize, 
    #                 percentPurityOfDS=percentPurityOfDS,
    #                 DL=DL, 
    #                 blenderSize=blenderSize, 
    #                 fillRatio=(fillRatio/100),
    #                 distribution= distribution)
    #       i += 1
    #     progbar.empty()
          
    #   flattenedResults = []
    #   for result in results:
    #     meanResults.append(np.mean(result))
    #     for val in result:
    #       flattenedResults.append(val)
  
      
    #   figure, ax = plt.subplots(figsize=(10,10))
    #   viz = ax.boxplot(meanResults)
    #   plt.title("Spread of mean assays for simulated blender with distribution = " + str(distribution))
    #   plt.ylabel("Mean Assay (%)")
    #   st.pyplot(figure)

    #   #show results table
    #   metrics = ["Min. Average Assay Observed (%)", 
    #              "Max. Average Assay Observed (%)", 
    #              "Standard Deviation", 
    #              "Min. Individual Assay Observed (%)", 
    #              "Max. Individual Assay Observed (%)"]
    #   values = [round(min(meanResults),2), 
    #             round(max(meanResults),2), 
    #             round(np.std(meanResults),2),
    #             round(min(flattenedResults),2),
    #             round(max(flattenedResults),2)]
      
    #   multipleRunResults = {"Metric": metrics,
    #                         "Value": values} 
    #   multipleRundf = pd.DataFrame(multipleRunResults)
    #   st.dataframe(multipleRundf)

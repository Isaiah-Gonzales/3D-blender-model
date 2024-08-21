import streamlit as st
from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib
import time
import numpy as np
import random
import math
import pandas as pd

def blender3D(blenderSize, fillRatio,thiefSize, distribution, DL=20, particleSize=100,percentPurityOfDS=100, clumpiness=0, clumpSize=1000):
  global placeholderaxes 
  global blender 
  global filledspace
  global top 
  global middle 
  global bottom 
  global TopSamplingArray
  global MidSamplingArray 
  global BotSamplingArray  
  
  #find out how many particles to simulate
  volPowder = blenderSize*fillRatio
  particlesizecm = particleSize/10000 #um to cm
  potentialParticles = blenderSize/particlesizecm
  particles = volPowder/particlesizecm
  particlesTop = potentialParticles *(5/6)
  placeholderaxes = int((particlesTop*(5/3)) ** (1/3))
  if placeholderaxes % 2 != 0:
      placeholderaxes+=1
    
  ##############################################################################################
  #create the shape which will represent the empty blender
  placeholderCube = np.zeros([placeholderaxes,placeholderaxes,placeholderaxes])
  slimmer = int(placeholderaxes/10)

  firstStepEnd = int(placeholderaxes*(4/10))
  secondStepEnd = int(placeholderaxes*(3/10))
  thirdStepEnd = int(placeholderaxes*(2/10))
  fourthStepEnd = int(placeholderaxes*(1/10))

  #create the blender shape out of the placeholder cube
  #top portion
  placeholderCube[:,:,firstStepEnd:placeholderaxes-1] = 1
  #second layer
  placeholderCube[slimmer:placeholderaxes-slimmer, slimmer:placeholderaxes-slimmer,secondStepEnd:firstStepEnd] = 1
  #third layer
  placeholderCube[slimmer*2:placeholderaxes-(slimmer*2), slimmer*2:placeholderaxes-(slimmer*2),thirdStepEnd:secondStepEnd] = 1
  #fourth layer
  placeholderCube[slimmer*3:placeholderaxes-(slimmer*3), slimmer*3:placeholderaxes-(slimmer*3),fourthStepEnd:thirdStepEnd] = 1
  #fifth layer
  placeholderCube[slimmer*4:placeholderaxes-(slimmer*4), slimmer*4:placeholderaxes-(slimmer*4),0:fourthStepEnd] = 1

  blender = placeholderCube.copy()

  ##############################################################################################
  #populate blender with powder
  filler = 0
  filledparticles = 0
  filledspace = blender.copy()
  while filledparticles < particles:
    filledparticles += np.sum(blender[:][:][filler])
    filler += 1

  #remove empty space on top of the filled space
  filledspace = filledspace[:,:,:filler+1]

  #portions to be populated with powder
  top = filledspace[:,:,firstStepEnd:]
  middle = filledspace[slimmer*2:placeholderaxes-(slimmer*2), slimmer*2:placeholderaxes-(slimmer*2),thirdStepEnd:secondStepEnd]
  bottom = filledspace[slimmer*4:placeholderaxes-(slimmer*4), slimmer*4:placeholderaxes-(slimmer*4),0:fourthStepEnd]

  ##############################################################################################
  #begin with all excipients
  top[:,:,:] = 0.00001 #close enough to zero where it won't affect simulation, but will still be displayed
  middle[:,:,:] = 0.00001
  bottom[:,:,:] = 0.00001

  numDSparticles = int(filledparticles * (DL/100))
  numDStop = int(np.prod(top.shape)*(DL/100))
  numDSmid = int(np.prod(middle.shape)*(DL/100))
  numDSbot = int(np.prod(bottom.shape)*(DL/100))
  portionTop = numDStop/numDSparticles
  portionMiddle = numDSmid/numDSparticles
  portionBottom = numDSbot/numDSparticles
  
  ##############################################################################################
  #now distribute DS
  if distribution == "unmixed":
    DSlevel = int(filler*(DL/100))
    if DSlevel > firstStepEnd:
      bottom[:,:,:] = percentPurityOfDS
      middle[:,:,:] = percentPurityOfDS
      top[:,:,:(DSlevel-firstStepEnd)+1] = percentPurityOfDS
    elif DSlevel > thirdStepEnd:
      bottom[:,:,:] = percentPurityOfDS
      if (DSlevel-secondStepEnd) < middle.shape[2]:
        middle[:,:,:(DSlevel-thirdStepEnd)+1] = percentPurityOfDS
      else:
        middle[:,:,:] = percentPurityOfDS
    else:
      if (DSlevel-fourthStepEnd) < bottom.shape[2]:
        bottom[:,:,:DSlevel+1] = percentPurityOfDS
      else:
        bottom[:,:,:] = percentPurityOfDS
        
  ##############################################################################################
  if distribution == "random":
    i = 0
    while i < numDStop:
      random_x = random.randint(0,top.shape[0]-1)
      random_z = random.randint(0,top.shape[1]-1)
      random_y = random.randint(0,top.shape[2]-1)
      if top[random_x,random_z,random_y] == 0.00001:
        top[random_x,random_z,random_y] = percentPurityOfDS
        i += 1
      else:
        pass
    i = 0
    while i < numDSmid:
      random_x = random.randint(0,middle.shape[0]-1)
      random_z = random.randint(0,middle.shape[1]-1)
      random_y = random.randint(0,middle.shape[2]-1)
      if middle[random_x,random_z,random_y] == 0.00001:
        middle[random_x,random_z,random_y] = percentPurityOfDS
        i += 1
      else:
        pass
    i = 0
    while i < numDSbot:
      random_x = random.randint(0,bottom.shape[0]-1)
      random_z = random.randint(0,bottom.shape[1]-1)
      random_y = random.randint(0,bottom.shape[2]-1)
      if bottom[random_x,random_z,random_y] == 0.00001:
        bottom[random_x,random_z,random_y] = percentPurityOfDS
        i += 1
      else:
        pass
    
  ##############################################################################################
  if distribution == "uniform":
    frequencyofDS = int(100/DL)
    flattenedtop = top.flatten()
    i=0
    while i < len(flattenedtop):
      if i%frequencyofDS == 0:
        flattenedtop[i] = percentPurityOfDS
        i += 1
      else:
        i+=1
    top = flattenedtop.reshape(top.shape)

    flattenedmiddle = middle.flatten()
    i=0
    while i < len(flattenedmiddle):
      if i%frequencyofDS == 0:
        flattenedmiddle[i] = percentPurityOfDS
        i += 1
      else:
        i+=1
    middle = flattenedmiddle.reshape(middle.shape)

    flattenedbottom = bottom.flatten()
    i=0
    while i < len(flattenedbottom):
      if i%frequencyofDS == 0:
        flattenedbottom[i] = percentPurityOfDS
        i += 1
      else:
        i+=1
    bottom = flattenedbottom.reshape(bottom.shape)
    
  ##############################################################################################
  if distribution == "poor":
    if clumpiness == 0:
      distribution = "random"
    else:
      clumpedParticles = int(numDSparticles*(clumpiness/10))
      numParticlesPerClump = int(clumpSize/particleSize)
      axisSizeclump = int(numParticlesPerClump**(1/3))
      if axisSizeclump > bottom.shape[2]:
        st.warning("Clump size too big in relation to blender for correct sampling, please reduce clump size or increase blender size.")
        return
      numClumps = int(clumpedParticles/numParticlesPerClump)
      #disperse clumps
      n = 0
      numClumpsTop = int(numClumps*portionTop)
      while n < numClumpsTop:
        random_x = random.randint(0,top.shape[0]-axisSizeclump)
        random_z = random.randint(0,top.shape[1]-axisSizeclump)
        random_y = random.randint(0,top.shape[2]-axisSizeclump)
        section = top[random_x:random_x+axisSizeclump,random_z:random_z+axisSizeclump,random_y:random_y+axisSizeclump]
        if any(percentPurityOfDS in sublist for sublist in section):
          pass
        else:
          top[random_x:random_x+axisSizeclump,random_z:random_z+axisSizeclump,random_y:random_y+axisSizeclump] = percentPurityOfDS
          n += 1
      n=0
      numClumpsMid = int(numClumps*portionMiddle)
      while n < numClumpsMid:
        random_x = random.randint(0,middle.shape[0]-axisSizeclump)
        random_z = random.randint(0,middle.shape[1]-axisSizeclump)
        random_y = random.randint(0,middle.shape[2]-axisSizeclump)
        section = middle[random_x:random_x+axisSizeclump,random_z:random_z+axisSizeclump,random_y:random_y+axisSizeclump]
        if any(percentPurityOfDS in sublist for sublist in section):
          pass
        else:
          middle[random_x:random_x+axisSizeclump,random_z:random_z+axisSizeclump,random_y:random_y+axisSizeclump] = percentPurityOfDS
          n += 1
      n=0
      numClumpsBot = int(numClumps*portionBottom)
      while n < numClumpsBot:
        random_x = random.randint(0,bottom.shape[0]-axisSizeclump)
        random_z = random.randint(0,bottom.shape[1]-axisSizeclump)
        random_y = random.randint(0,bottom.shape[2]-axisSizeclump)
        section = bottom[random_x:random_x+axisSizeclump,random_z:random_z+axisSizeclump,random_y:random_y+axisSizeclump]
        if any(percentPurityOfDS in sublist for sublist in section):
          pass
        else:
          bottom[random_x:random_x+axisSizeclump,random_z:random_z+axisSizeclump,random_y:random_y+axisSizeclump] = percentPurityOfDS
          n += 1

      remainingParticles = (numDStop+numDSmid+numDSbot) - clumpedParticles
      remainderTop = int(remainingParticles*portionTop)
      remainderMid = int(remainingParticles*portionMiddle)
      reminderBot = int(remainingParticles*portionBottom)

      i=0
      while i < remainderTop:
        random_x = random.randint(0,top.shape[0]-1)
        random_z = random.randint(0,top.shape[1]-1)
        random_y = random.randint(0,top.shape[2]-1)
        if top[random_x,random_z,random_y] == 0.00001:
          top[random_x,random_z,random_y] = percentPurityOfDS
          i += 1
        else:
          pass
      i=0
      while i < remainderMid:
        random_x = random.randint(0,middle.shape[0]-1)
        random_z = random.randint(0,middle.shape[1]-1)
        random_y = random.randint(0,middle.shape[2]-1)
        if middle[random_x,random_z,random_y] == 0.00001:
          middle[random_x,random_z,random_y] = percentPurityOfDS
          i += 1
        else:
          pass
      i=0
      while i < reminderBot:
        random_x = random.randint(0,bottom.shape[0]-1)
        random_z = random.randint(0,bottom.shape[1]-1)
        random_y = random.randint(0,bottom.shape[2]-1)
        if bottom[random_x,random_z,random_y] == 0.00001:
          bottom[random_x,random_z,random_y] = percentPurityOfDS
          i += 1
        else:
          pass
      
  ##############################################################################################
  #sampling
  results = []
  toSample = {}
  numSamples = 4
  currentSample = 0
  thiefDimensions = int((thiefSize/particlesizecm)**(1/3))
  if thiefDimensions > bottom.shape[2]:
    st.warning("thief size too large in comparison to blender. Please increase blender size or decrease thief size")
    return
  
  topfirststart = slimmer*4
  topfirstend = (slimmer*4) + thiefDimensions
  topsecondstart = placeholderaxes-(thiefDimensions+(slimmer*4))
  topsecondend = placeholderaxes-(slimmer*4)

  midfirststart = slimmer*2
  midfirstend = (slimmer*2) + thiefDimensions
  midsecondstart = middle.shape[1]-(thiefDimensions+(slimmer*2))
  midsecondend = middle.shape[1]-(slimmer*2)

  botfirststart = 0
  botfirstend = thiefDimensions
  botsecondstart = bottom.shape[1]-(thiefDimensions)
  botsecondend = bottom.shape[1]

  #top sampling
  TopSamplingArray = np.zeros(top.shape)
  #top-1
  TopSamplingArray[topfirststart:topfirstend,topfirststart:topfirstend,top.shape[2]-thiefDimensions:top.shape[2]] = 1
  results.append(np.mean(top[topfirststart:topfirstend,topfirststart:topfirstend,top.shape[2]-thiefDimensions:top.shape[2]]))
  #top-2
  TopSamplingArray[topsecondstart:topsecondend,topfirststart:topfirstend,top.shape[2]-thiefDimensions:top.shape[2]] = 1
  results.append(np.mean(top[topsecondstart:topsecondend,topfirststart:topfirstend,top.shape[2]-thiefDimensions:top.shape[2]]))
  #top-3
  TopSamplingArray[topfirststart:topfirstend,topsecondstart:topsecondend,top.shape[2]-thiefDimensions:top.shape[2]] = 1
  results.append(np.mean(top[topfirststart:topfirstend,topsecondstart:topsecondend,top.shape[2]-thiefDimensions:top.shape[2]]))
  #top-4
  TopSamplingArray[topsecondstart:topsecondend,topsecondstart:topsecondend,top.shape[2]-thiefDimensions:top.shape[2]] = 1
  results.append(np.mean(top[topsecondstart:topsecondend,topsecondstart:topsecondend,top.shape[2]-thiefDimensions:top.shape[2]]))
  
  #middle sampling
  MidSamplingArray = np.zeros(middle.shape)
  #mid-1
  MidSamplingArray[midfirststart:midfirstend,midfirststart:midfirstend,middle.shape[2]-thiefDimensions:middle.shape[2]] = 1
  results.append(np.mean(middle[midfirststart:midfirstend,midfirststart:midfirstend,middle.shape[2]-thiefDimensions:middle.shape[2]]))
  #mid-2
  MidSamplingArray[midsecondstart:midsecondend,midfirststart:midfirstend,middle.shape[2]-thiefDimensions:middle.shape[2]] = 1
  results.append(np.mean(middle[midsecondstart:midsecondend,midfirststart:midfirstend,middle.shape[2]-thiefDimensions:middle.shape[2]]))
  #mid-3
  MidSamplingArray[midfirststart:midfirstend,midsecondstart:midsecondend,middle.shape[2]-thiefDimensions:middle.shape[2]] = 1
  results.append(np.mean(middle[midfirststart:midfirstend,midsecondstart:midsecondend,middle.shape[2]-thiefDimensions:middle.shape[2]]))
  #mid-4
  MidSamplingArray[midsecondstart:midsecondend,midsecondstart:midsecondend,middle.shape[2]-thiefDimensions:middle.shape[2]] = 1
  results.append(np.mean(middle[midsecondstart:midsecondend,midsecondstart:midsecondend,middle.shape[2]-thiefDimensions:middle.shape[2]]))

  #bottom sampling
  BotSamplingArray = np.zeros(bottom.shape)
  #bot-1
  BotSamplingArray[botfirststart:botfirstend,botfirststart:botfirstend,bottom.shape[2]-thiefDimensions:bottom.shape[2]] = 1
  results.append(np.mean(bottom[botfirststart:botfirstend,botfirststart:botfirstend,bottom.shape[2]-thiefDimensions:bottom.shape[2]]))
  #bot-2
  BotSamplingArray[botsecondstart:botsecondend,botfirststart:botfirstend,bottom.shape[2]-thiefDimensions:bottom.shape[2]] = 1
  results.append(np.mean(bottom[botsecondstart:botsecondend,botfirststart:botfirstend,bottom.shape[2]-thiefDimensions:bottom.shape[2]]))
  #bot-3
  BotSamplingArray[botfirststart:botfirstend,botsecondstart:botsecondend,bottom.shape[2]-thiefDimensions:bottom.shape[2]] = 1
  results.append(np.mean(bottom[botfirststart:botfirstend,botsecondstart:botsecondend,bottom.shape[2]-thiefDimensions:bottom.shape[2]]))
  #bot-4
  BotSamplingArray[botsecondstart:botsecondend,botsecondstart:botsecondend,bottom.shape[2]-thiefDimensions:bottom.shape[2]] = 1
  results.append(np.mean(bottom[botsecondstart:botsecondend,botsecondstart:botsecondend,bottom.shape[2]-thiefDimensions:bottom.shape[2]]))

  #resize top, middle, and bottom so they'll be correctly displayed      
  top = np.pad(top,((0,0),(0,0),(firstStepEnd,0)), mode='constant')
  middle = np.pad(middle,((slimmer*2,slimmer*2),(slimmer*2,slimmer*2),(thirdStepEnd,0)))
  bottom = np.pad(bottom,((slimmer*4,slimmer*4),(slimmer*4,0),(0,0)))

  TopSamplingArray = np.pad(TopSamplingArray,((0,0),(0,0),(firstStepEnd,0)), mode='constant')
  MidSamplingArray = np.pad(MidSamplingArray,((slimmer*2,slimmer*2),(slimmer*2,slimmer*2),(thirdStepEnd,0)))
  BotSamplingArray = np.pad(BotSamplingArray,((slimmer*4,slimmer*4),(slimmer*4,0),(0,0)))

  #adjust results so they're in percent of expected assay
  percentAssays = []
  for results in results:
    assay = ((results/DL)*100)
    percentAssays.append(round(assay,2))
  
  #display results table
  positions = ["top-1", "top-2", "top-3", "top-4",
               "mid-1", "mid-2", "mid-3", "mid-4",
               "bot-1", "bot-2", "bot-3", "bot-4"]
  data1 = {'positions': positions,
           'assay (%)': percentAssays}
  df = pd.DataFrame(data1)
  st.dataframe(data=df)
  
def displayBlender(placeholderaxes, blender,filledspace, top, middle, bottom, TopSamplingArray, MidSamplingArray, BotSamplingArray):
  progbar = st.progress(0.0, text = "Initializing...")

  ##############################################################################################
  fig = plt.figure(figsize=(33,10))
  ax = fig.add_subplot(1, 3, 1, projection='3d')
  ax.set_title("Blender shape vs. volume filled with powder")
  ax.voxels(blender, facecolors='lightgray', alpha=0.5)
  ax.voxels(filledspace, facecolors='white')
  ax.set_xlim(0,placeholderaxes+2)
  ax.set_ylim(0,placeholderaxes+2)
  ax.set_zlim(0,placeholderaxes+2)
  ax.view_init(elev=10)
  #axis labels
  dimensionsinmL = round(placeholderaxes*(particleSize/1000),2)
  #ax.set_xticks([])
  ax.set_xlabel("tick marks show particle count equivalent to: " + str(dimensionsinmL) + " cm")
  #ax.set_yticks([])
  ax.set_ylabel(str(dimensionsinmL) + " cm")
  #ax.set_zticks([])
  ax.set_zlabel(str(dimensionsinmL) + " cm")
  progbar.progress(0.25, text = "Plot 1 of 3 generated")
  
  ##############################################################################################
  ax2 = fig.add_subplot(1, 3, 2, projection='3d')
  ax2.set_title("Top, middle, and bottom sample locations with DS distribution = " + str(distribution))
  cmap = plt.cm.binary.copy()
  topcolors = cmap(top)
  ax2.voxels(top, facecolors=topcolors)
  midcolors = cmap(middle)
  ax2.voxels(middle, facecolors=midcolors)
  botcolors = cmap(bottom)
  ax2.voxels(bottom, facecolors=botcolors)
  ax2.voxels(TopSamplingArray, facecolors='b')
  ax2.voxels(MidSamplingArray, facecolors='y')
  ax2.voxels(BotSamplingArray, facecolors='r')
  ax2.set_xlim(0,placeholderaxes+2)
  ax2.set_ylim(0,placeholderaxes+2)
  ax2.set_zlim(0,placeholderaxes+2)
  ax2.view_init(elev=10)
  ax2.set_xlabel("tick marks show particle count equivalent to: " + str(dimensionsinmL) + " cm")
  ax2.set_ylabel(str(dimensionsinmL) + " cm")
  ax2.set_zlabel(str(dimensionsinmL) + " cm")
  progbar.progress(0.50, text = "Plot 2 of 3 generated")
  
  ##############################################################################################
  ax3 = fig.add_subplot(1, 3, 3, projection='3d')
  ax3.set_title("Samples extracted")
  ax3.voxels(top,alpha = 0.15)
  ax3.voxels(middle,alpha = 0.15)
  ax3.voxels(bottom,alpha = 0.15)
  ax3.voxels(TopSamplingArray, facecolors="b")
  ax3.voxels(MidSamplingArray, facecolors="y")
  ax3.voxels(BotSamplingArray, facecolors="r")
  ax3.set_xlim(0,placeholderaxes+2)
  ax3.set_ylim(0,placeholderaxes+2)
  ax3.set_zlim(0,placeholderaxes+2)
  ax3.view_init(elev=10)
  ax3.set_xlabel("tick marks show particle count equivalent to: " + str(dimensionsinmL) + " cm")
  ax3.set_ylabel(str(dimensionsinmL) + " cm")
  ax3.set_zlabel(str(dimensionsinmL) + " cm")    
  norm = matplotlib.colors.Normalize(vmin=0, vmax=percentPurityOfDS)
  m=plt.cm.ScalarMappable(cmap=plt.cm.binary, norm=norm)
  m.set_array([])
  plt.colorbar(m, label='DS %')
  progbar.progress(0.75, text = "Plot 3 of 3 generated")

  ##############################################################################################
  plt.tight_layout()
  st.pyplot(fig)
  progbar.empty()

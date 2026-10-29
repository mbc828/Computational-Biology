# Public source snapshot: replace path placeholders before running.
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(r'<PROJECT_ROOT>\outputs\manuscript_v21_github_release\figure_consolidation_v82')
MEDIA=ROOT/'work/unpacked/word/media'
OUT=ROOT/'composite_figures'
OUT.mkdir(exist_ok=True)
FONT=str(Path(r'<LOCAL_DRIVE>/Windows\Fonts\arialbd.ttf'))

def im(n):
    return Image.open(MEDIA/f'image{n}.png').convert('RGB')

def relabel(img, replacements, size=70):
    """ replacements: (x,y,old_box_w,old_box_h,new_letter) in current pixels """
    d=ImageDraw.Draw(img)
    f=ImageFont.truetype(FONT,size)
    for x,y,w,h,letter in replacements:
        d.rectangle((x-8,y-5,x+w+8,y+h+5),fill='white')
        d.text((x,y),letter,font=f,fill='#20252b')
    return img

def add_label(img,letter,x=25,y=20,size=75):
    return relabel(img,[(x,y,size,size,letter)],size)

def canvas(rows,width=3600,margin=40,gap=45):
    """Each row is a list of images. Images are scaled to equal cell widths."""
    placed=[]
    y=margin
    for images in rows:
        cell_w=(width-2*margin-gap*(len(images)-1))//len(images)
        scaled=[]
        for image in images:
            h=round(image.height*cell_w/image.width)
            scaled.append(image.resize((cell_w,h),Image.Resampling.LANCZOS))
        row_h=max(x.height for x in scaled)
        placed.append((y,cell_w,scaled,row_h))
        y+=row_h+gap
    sheet=Image.new('RGB',(width,y-gap+margin),'white')
    for y0,cell_w,images,row_h in placed:
        for j,image in enumerate(images):
            x=margin+j*(cell_w+gap)+(cell_w-image.width)//2
            sheet.paste(image,(x,y0))
    return sheet

def save(sheet,name):
    path=OUT/f'{name}.png'
    sheet.save(path,optimize=True,dpi=(300,300))
    prev=sheet.copy()
    prev.thumbnail((1800,1800),Image.Resampling.LANCZOS)
    prev.save(OUT/f'{name}_preview.jpg',quality=90,optimize=True)
    print(name,sheet.size,path.stat().st_size,'preview',prev.size)

def with_header(img,text,height=150,size=62):
    panel=Image.new('RGB',(img.width,img.height+height),'white')
    panel.paste(img,(0,height))
    d=ImageDraw.Draw(panel)
    d.text((35,25),text,font=ImageFont.truetype(FONT,size),fill='#20252b')
    return panel

# Figure 1: corrected patient transcriptome, A-F.
a9=im(9).crop((0,330,4251,2500))
go=im(10)
# Keep the full common colour scale with the UP panel instead of splitting it.
cplot=go.crop((0,620,3700,2600))
cb=go.crop((3700,950,4251,4050))
cb=cb.resize((round(cb.width*cplot.height/cb.height),cplot.height),Image.Resampling.LANCZOS)
cgo=Image.new('RGB',(cplot.width+cb.width,cplot.height),'white')
cgo.paste(cplot,(0,0)); cgo.paste(cb,(cplot.width,0))
c10=with_header(cgo,'C  UP: top 10 of 664 significant terms',140,66)
dplot=go.crop((0,2850,3700,4600))
dgo=Image.new('RGB',(cgo.width,dplot.height),'white')
dgo.paste(dplot,(0,0))
d10=with_header(dgo,'D  DOWN: all 7 significant terms',140,66)
e11=with_header(im(11).crop((160,650,2140,2200)),'E  TLR2 / TLR4 / MyD88',140,68)
f11=with_header(im(11).crop((160,2400,2140,3800)),'F  Canonical NF-κB',140,68)
save(canvas([[a9],[c10,d10],[e11,f11]]),'Figure_1_clinical_transcriptome')

# Figure 2: workflow and candidate nomination, A-D.
a1=add_label(im(1),'A',30,25,80)
b6=im(6).crop((3975,0,7775,3700))
bd=ImageDraw.Draw(b6)
bd.rectangle((0,3100,1100,3600),fill='white')
bd.ellipse((70,3340,140,3410),fill='#dceaf5',outline='#6b8ba4',width=3)
bd.text((175,3342),'Scutellaria compound',font=ImageFont.truetype(str(Path(r'<LOCAL_DRIVE>/Windows\Fonts\arial.ttf')),52),fill='#20252b')
# Redraw the compound centrality panel from the values displayed in original Figure 6C.
compounds=['Berberine','Baicalein','Emodin','Baicalin','Wogonin','Apigenin','Coptisine','Palmatine','Oroxylin A','Jatrorrhizine','Wogonoside','Chrysin','Rhein','Aloe-emodin','Chrysophanol','Magnoflorine']
centralities=[0.490,0.480,0.440,0.410,0.390,0.370,0.350,0.340,0.330,0.320,0.310,0.300,0.290,0.280,0.270,0.250]
herbs=['Coptis','Scutellaria','Rheum','Scutellaria','Scutellaria','Scutellaria','Coptis','Coptis','Scutellaria','Coptis','Scutellaria','Scutellaria','Rheum','Rheum','Rheum','Coptis']
colors={'Scutellaria':'#2878a5','Coptis':'#55b89b','Rheum':'#c28b16'}
figc,axc=plt.subplots(figsize=(6.7,6.0),dpi=300)
y=np.arange(len(compounds))[::-1]
axc.hlines(y,0,centralities,color='#bbbbbb',lw=0.7)
for yi,v,h in zip(y,centralities,herbs):
    axc.plot(v,yi,'o',color=colors[h],ms=5)
    axc.text(v+0.012,yi,f'{v:.3f}',va='center',fontsize=7.5,color=colors[h])
axc.set_yticks(y); axc.set_yticklabels(compounds,fontsize=8)
axc.set_xlim(0,0.56); axc.set_xlabel('Normalized centrality score',fontsize=9)
axc.tick_params(axis='x',labelsize=8)
for s in ('top','right'): axc.spines[s].set_visible(False)
axc.set_title('C  Compound centrality and connectivity',loc='left',fontweight='bold',fontsize=11)
figc.tight_layout()
comp_path=OUT/'compound_centrality.png'
figc.savefig(comp_path,bbox_inches='tight',pad_inches=0.08,facecolor='white')
plt.close(figc)
c6=Image.open(comp_path).convert('RGB')
d7plot=im(7).crop((1400,150,3602,1658))
ImageDraw.Draw(d7plot).rectangle((0,450,25,700),fill='white')
d7=with_header(d7plot,'D  Top target centrality scores',130,60)
save(canvas([[a1,b6],[c6,d7]]),'Figure_2_candidate_nomination')

# Figure 3: docking and IKKbeta MD, A-F.
a13=im(13)
# Extract the IKKβ energy panel and restore its shared x-axis below the legend.
ikk=im(14).crop((0,700,4251,2100))
ImageDraw.Draw(ikk).rectangle((450,1270,1300,1400),fill='white')
xaxis=im(14).crop((500,3520,4251,3800))
xd=ImageDraw.Draw(xaxis)
xd.line((200,0,3500,0),fill='#20252b',width=3)
for xx in (460,1080,1670,2280,2890,3500):
    xd.line((xx,0,xx,20),fill='#20252b',width=3)
ikk_with_axis=Image.new('RGB',(4251,ikk.height+xaxis.height),'white')
ikk_with_axis.paste(ikk,(0,0)); ikk_with_axis.paste(xaxis,(500,ikk.height))
c14=with_header(ikk_with_axis,'C  IKKβ',140,72)
# Preserve the original three-panel row and change only its panel letters.
def16=relabel(im(16),[(195,0,30,55,'D'),(1205,0,30,55,'E'),(2347,0,30,55,'F')],40)
save(canvas([[a13],[c14],[def16]]),'Figure_3_IKKb_MD_energy_fingerprint')

# Figure 4F: explicit mutant MD means transcribed from Table A7-1.
plt.rcParams.update({'font.family':'Arial','font.size':9})
fig,axes=plt.subplots(1,2,figsize=(11.8,3.5),dpi=300,gridspec_kw={'width_ratios':[2.3,1]})
sets=[('I165A',['Baicalein','Oroxylin A','Berberine','Wogonin','Apigenin','Coptisine'],[6.46,5.98,7.36,0.61,0.07,3.46]),('D166A',['Baicalein','Apigenin'],[-4.94,-3.28])]
for ax,(title,names,vals) in zip(axes,sets):
    x=np.arange(len(names))
    colors=['#2878a5' if v>=0 else '#c96a3b' for v in vals]
    bars=ax.bar(x,vals,color=colors,width=0.66)
    ax.axhline(0,color='#333333',lw=0.8)
    ax.set_title(title,loc='left',fontweight='bold',fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(names,rotation=32 if len(names)>2 else 0,ha='right' if len(names)>2 else 'center',fontsize=8)
    ax.set_ylabel(r'$\Delta\Delta G$ (kcal/mol)')
    ax.grid(axis='y',color='#dddddd',lw=0.5,alpha=0.7)
    ax.set_axisbelow(True)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    for b,v in zip(bars,vals):
        ax.text(b.get_x()+b.get_width()/2,v+(0.25 if v>=0 else -0.35),f'{v:+.2f}',ha='center',va='bottom' if v>=0 else 'top',fontsize=8)
    ax.set_ylim(-6.4,9.0)
fig.text(0.995,0.01,'Explicit mutant MD; n = 3 trajectories per system',ha='right',fontsize=8,color='#444444')
fig.tight_layout(rect=(0,0.04,1,1),w_pad=2.0)
mut=OUT/'mutation_panel.png'
fig.savefig(mut,bbox_inches='tight',pad_inches=0.08,facecolor='white')
plt.close(fig)
fmut=add_label(Image.open(mut).convert('RGB'),'F',15,10,65)

a19=im(19)
c20=relabel(im(20),[(0,0,95,105,'C'),(1450,0,95,105,'D'),(2890,0,95,105,'E')],82)
save(canvas([[a19],[c20],[fmut]]),'Figure_4_binding_pose_mutation')
print('done')

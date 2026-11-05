"""Reproduce manuscript figures. Run: python analysis.py [input.xlsx] [output_dir]."""
from pathlib import Path
import sys, json, platform, os, argparse
os.environ.setdefault('MPLCONFIGDIR','/tmp/rasc_matplotlib')
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docx import Document
from docx.shared import Inches

parser = argparse.ArgumentParser(description='Generate manuscript CUSUM figures from a local 48-case Excel dataset.')
parser.add_argument('input_xlsx', type=Path, help='Local workbook containing Sheet1; never commit patient data.')
parser.add_argument('output_dir', type=Path, nargs='?', default=Path('outputs'), help='Output directory (default: outputs)')
args = parser.parse_args()
SOURCE = args.input_xlsx
OUT = args.output_dir
OUT.mkdir(parents=True, exist_ok=True)
d = pd.read_excel(SOURCE, sheet_name='Sheet1')
assert len(d)==48 and d['病例编号'].is_unique
assert d['手术日期'].is_monotonic_increasing
assert d['术者_文件现有标识'].eq('术者1').all()
assert d['手术入路'].eq('机器人辅助').all()
n=len(d); case=np.arange(1,n+1)
cols=['BMI_kg_m2','同期子宫切除','后壁或会阴修补','抗尿失禁手术']
assert not d[cols+['总手术时间_min']].isna().any().any()
assert np.array_equal(d['后壁或会阴修补'],np.maximum(d['后壁修补'],d['会阴修补或成形']))
y=d['总手术时间_min'].to_numpy(float)
X=np.c_[np.ones(n),d[cols].to_numpy(float)]

def fit(x, z):
    assert np.linalg.matrix_rank(x)==x.shape[1]
    b=np.linalg.lstsq(x,z,rcond=None)[0]
    e=z-x@b
    inv=np.linalg.inv(x.T@x)
    h=np.sum((x@inv)*x,axis=1)
    cov=inv@(x.T@((e/(1-h))[:,None]**2*x))@inv
    se=np.sqrt(np.diag(cov)); df=len(z)-x.shape[1]
    # t-based HC3 intervals; normal-based values also exported to expose convention.
    ci=stats.t.ppf(.975,df)*se
    r2=1-np.sum(e**2)/np.sum((z-z.mean())**2)
    return dict(b=b,e=e,pred=x@b,se=se,lo=b-ci,hi=b+ci,
                p=2*stats.t.sf(abs(b/se),df),p_normal=2*stats.norm.sf(abs(b/se)),
                r2=r2,adj_r2=1-(1-r2)*(len(z)-1)/df,
                loo_mae=np.mean(abs(e/(1-h))))

m=fit(X,y); residual=m['e']; raw=y-y.mean()
assert abs(residual.sum())<1e-7
cr=np.cumsum(raw); ca=np.cumsum(residual)
assert np.isclose(cr.max(),497.875) and cr.argmax()+1==31
assert np.isclose(ca.max(),488.3710117886646) and ca.argmax()+1==21

def permute(v, seed):
    rng=np.random.default_rng(seed)
    curves=np.cumsum(np.array([rng.permutation(v) for _ in range(10000)]),axis=1)
    observed=np.max(abs(np.cumsum(v)))
    exceed=int(np.sum(np.max(abs(curves),axis=1)>=observed))
    return np.quantile(curves,[.025,.975],axis=0),(exceed+1)/10001,exceed

band_r,pr,er=permute(raw,20261105)
band_a,pa,ea=permute(residual,20261106)

# POP-Q stage from maximal measured leading edge, relative to TVL.
# Two structured fields are missing. Original strings show B+3 (case 18,
# ambiguous point label) and a+3 (case 36). Available points already give
# stages III and IV; these original +3 values cannot change those stages.
lead=d[['Aa_cm','Ba_cm','C_cm','Ap_cm','Bp_cm']].max(axis=1)
stage=np.select([lead < -1,lead<=1,lead<d['TVL_cm']-2],[1,2,3],default=4)
assert pd.Series(stage).value_counts().to_dict()=={3:30,4:9,2:9}
Xp=np.c_[X,stage]; mp=fit(Xp,y); cp=np.cumsum(mp['e'])

loo=[]; loo_curves=[]
for i in range(n):
    keep=case!=i+1; mm=fit(X[keep],y[keep]); cc=np.cumsum(mm['e'])
    j=int(np.argmax(cc)); orig=int(case[keep][j])
    loo.append(dict(omitted_case=i+1,peak_original_case=orig,peak_reduced_index=j+1,peak_min=float(cc[j])))
    cc_full=np.full(n,np.nan);cc_full[keep]=cc;loo_curves.append(cc_full)
loo=pd.DataFrame(loo)

smooth={}; alternatives=[]
for label,v in [('Traditional',y),('Adjusted',residual)]:
    for edge in ['complete','shortened']:
        sm=pd.Series(v).rolling(5,center=True,min_periods=5 if edge=='complete' else 1).mean()
        cc=(sm-sm.mean()).cumsum()
        alternatives.append(dict(curve=label,edge=edge,recentered=True,peak_original_case=int(cc.idxmax()+1),peak_min=float(cc.max())))
        if edge=='complete': smooth[label]=cc.to_numpy()

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
 'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'axes.titlesize':11,'axes.labelsize':10,'legend.fontsize':8,'savefig.facecolor':'white'})
BLUE='#236B8E'; ORANGE='#C65D30'; GREEN='#32816B'; GRAY='#B7C7D4'
def base(ax, ylabel='Cumulative deviation (min)',zero=True):
    ax.set(xlim=(0,49),xlabel='Consecutive case number',ylabel=ylabel)
    ax.set_xticks([0,10,20,30,40,48]);ax.grid(axis='y',alpha=.17)
    if zero: ax.axhline(0,color='.45',lw=.8)
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=400,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.tiff'),dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)
def peak(ax,c,color):
    j=int(np.nanargmax(c));ax.scatter(j+1,c[j],color=color,s=32,zorder=4)
    ax.annotate(f'Case {j+1}: {c[j]:.2f} min',xy=(j+1,c[j]),xytext=(0,14),
                textcoords='offset points',ha='center',fontsize=9,color=color)

fig,ax=plt.subplots(figsize=(8,4.7),layout='constrained')
ax.plot(case,y,'o-',color=BLUE,ms=3.5,lw=.9,label='Observed total operative time')
ax.plot(case,m['pred'],color=ORANGE,lw=1.6,ls='--',label='Case-mix predicted time (no case number)')
ax.axhline(y.mean(),color='.4',lw=1,ls=':',label=f'Cohort mean: {y.mean():.2f} min')
base(ax,'Total operative time (min)',False);ax.legend(loc='upper right',frameon=False)
save(fig,'Figure_1_operative_time')

fig,axs=plt.subplots(1,2,figsize=(11,4.8),sharey=True,layout='constrained')
for ax,c,band,p,color,title in [(axs[0],cr,band_r,pr,BLUE,'A  Traditional CUSUM'),(axs[1],ca,band_a,pa,ORANGE,'B  Case-mix adjusted CUSUM')]:
    ax.fill_between(case,band[0],band[1],color=GRAY,alpha=.55,label='Pointwise 95% permutation reference interval')
    ax.plot(np.r_[0,case],np.r_[0,c],color=color,lw=2,label='Observed CUSUM')
    peak(ax,c,color);base(ax);ax.set_title(title,loc='left',pad=12)
    ax.text(.04,.07,f'Global permutation P = {p:.3f}\n10,000 permutations',transform=ax.transAxes,fontsize=9)
    ax.set_ylim(min(band_r.min(),band_a.min())-35,max(cr.max(),ca.max())+100)
axs[0].legend(loc='lower left',bbox_to_anchor=(0,.15),frameon=False,fontsize=7)
save(fig,'Figure_2_CUSUM_permutation')

fig,axs=plt.subplots(2,2,figsize=(11,8.1),layout='constrained')
ax=axs[0,0]
for c in loo_curves: ax.plot(case,c,color=GRAY,lw=.8,alpha=.7)
ax.plot(case,ca,color=ORANGE,lw=2,label='Full 48-case model');base(ax)
ax.set_title('A  Leave-one-out adjusted CUSUM',loc='left');ax.legend(frameon=False)
ax=axs[0,1];counts=loo.peak_original_case.value_counts().sort_index()
ax.bar(counts.index.astype(str),counts.values,color=ORANGE,width=.6)
for i,v in enumerate(counts.values):ax.text(i,v+.5,str(v),ha='center')
ax.set(xlabel='Peak location (original case number)',ylabel='Number of leave-one-out refits',ylim=(0,52))
ax.set_title('B  Peak-location stability (48 refits)',loc='left');ax.grid(axis='y',alpha=.17)
ax=axs[1,0];ax.plot(case,ca,color=ORANGE,lw=2,label='Primary adjustment')
ax.plot(case,cp,color=GREEN,lw=1.7,ls='--',label='Primary adjustment + POP-Q stage')
j=int(cp.argmax());ax.scatter(j+1,cp[j],color=GREEN,s=32,zorder=4)
ax.annotate(f'Case {j+1}: {cp[j]:.2f} min',xy=(j+1,cp[j]),xytext=(.05,.92),textcoords='axes fraction',color=GREEN,fontsize=9,arrowprops={'arrowstyle':'-','color':GREEN,'lw':.8})
base(ax);ax.set_title('C  Additional adjustment for POP-Q stage',loc='left');ax.legend(frameon=False,loc='lower left');ax.margins(y=.2)
ax=axs[1,1]
for label,color in [('Traditional',BLUE),('Adjusted',ORANGE)]:
    cc=smooth[label];j=int(np.nanargmax(cc))
    ax.plot(case,cc,color=color,lw=2,label=label);ax.scatter(j+1,cc[j],color=color,s=32,zorder=4)
    ax.text(.04 if label=='Adjusted' else .96,.94,f'Case {j+1}: {cc[j]:.2f} min',transform=ax.transAxes,ha='left' if label=='Adjusted' else 'right',color=color,fontsize=8)
base(ax);ax.set_title('D  Centered 5-case moving-average CUSUM',loc='left');ax.legend(frameon=False,loc='lower left');ax.margins(y=.25)
save(fig,'Figure_S1_sensitivity')

trends=[]
for name,xx in [('Unadjusted',np.c_[np.ones(n),case]),('Primary adjusted',np.c_[X,case]),('POP-Q adjusted',np.c_[Xp,case])]:
    mm=fit(xx,y);trends.append(dict(model=name,beta_per_case=mm['b'][-1],hc3_se=mm['se'][-1],ci_low_t=mm['lo'][-1],ci_high_t=mm['hi'][-1],p_t=mm['p'][-1],ci_low_normal=mm['b'][-1]-stats.norm.ppf(.975)*mm['se'][-1],ci_high_normal=mm['b'][-1]+stats.norm.ppf(.975)*mm['se'][-1],p_normal=mm['p_normal'][-1]))
trends=pd.DataFrame(trends)
curves=pd.DataFrame({'case_number':case,'operative_time_min':y,'predicted_time_min':m['pred'],'residual_min':residual,
 'traditional_CUSUM':cr,'adjusted_CUSUM':ca,'traditional_perm_lower':band_r[0],'traditional_perm_upper':band_r[1],
 'adjusted_perm_lower':band_a[0],'adjusted_perm_upper':band_a[1],'POP_Q_stage':stage,'POP_adjusted_CUSUM':cp,
 'smoothed_traditional_CUSUM':smooth['Traditional'],'smoothed_adjusted_CUSUM':smooth['Adjusted']})
coef=pd.DataFrame({'term':['Intercept','BMI','Hysterectomy','Posterior/perineal repair','Anti-incontinence procedure'],
 'coefficient':m['b'],'HC3_SE':m['se'],'HC3_CI_low_t':m['lo'],'HC3_CI_high_t':m['hi'],'HC3_P_t':m['p']})
diff=d['总手术时间_min'].ne(d['手术时长_原值'])
audit=pd.DataFrame({'case_number':case,'analysis_time_min':y,'original_time_min':d['手术时长_原值'],
 'time_difference_flag':diff,'same_date_order_requires_confirmation':d['手术日期'].duplicated(keep=False),
 'POP_structured_point_missing':d[['Aa_cm','Ba_cm','C_cm','Ap_cm','Bp_cm']].isna().any(axis=1)})
with pd.ExcelWriter(OUT/'analysis_results.xlsx',engine='openpyxl') as w:
    curves.to_excel(w,sheet_name='Figure_data',index=False)
    coef.to_excel(w,sheet_name='Case_mix_model',index=False)
    trends.to_excel(w,sheet_name='Time_trend_HC3',index=False)
    loo.to_excel(w,sheet_name='Leave_one_out',index=False)
    pd.DataFrame(alternatives).to_excel(w,sheet_name='Smoothing_definitions',index=False)
    audit.to_excel(w,sheet_name='Data_audit',index=False)
    pd.DataFrame([dict(curve='Traditional',seed=20261105,permutations=10000,exceedances=er,p=pr),dict(curve='Adjusted',seed=20261106,permutations=10000,exceedances=ea,p=pa)]).to_excel(w,sheet_name='Permutation_tests',index=False)

summary={'n':n,'mean_time':y.mean(),'sd_time':y.std(ddof=1),'raw_peak_case':int(cr.argmax()+1),'raw_peak':cr.max(),
 'adjusted_peak_case':int(ca.argmax()+1),'adjusted_peak':ca.max(),'permutation_p_raw':pr,'permutation_p_adjusted':pa,
 'POP_peak_case':int(cp.argmax()+1),'POP_peak':cp.max(),'loo_peak_original_counts':{str(k):int(v) for k,v in counts.items()},
 'model_r2':m['r2'],'model_adjusted_r2':m['adj_r2'],'model_loo_mae':m['loo_mae'],
 'trend':trends.to_dict('records'),'smoothing':alternatives,
 'versions':dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__)}
(OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')

captions=[
('Figure 1. Observed and case-mix predicted total operative time.',
 'Observed operative time and fitted values from an ordinary least-squares model including an intercept, BMI, concomitant hysterectomy, posterior vaginal wall and/or perineal repair, and anti-incontinence surgery. Case number was excluded from the prediction model. Cases follow the existing chronological order in Sheet1; within-date order was retained.'),
('Figure 2. Traditional and case-mix adjusted CUSUM curves.',
 f'(A) Traditional CUSUM of operative time minus the 48-case mean. (B) Cumulative residuals from the case-mix model. Shading indicates pointwise 95% permutation reference intervals from 10,000 permutations, not simultaneous confidence bands or proficiency thresholds. Two-sided global tests use the maximum absolute cumulative deviation; Monte Carlo P = (b+1)/(10000+1). Seeds: 20261105 and 20261106. Peaks occurred at cases 31 (497.88 min) and 21 (488.37 min), respectively. Global P values were {pr:.3f} and {pa:.3f}. Residual permutation treats fitted residuals as exchangeable and does not repeat model fitting; results are exploratory under this assumption. Peak locations describe cumulative deviations and do not establish proficiency.'),
('Supplementary Figure S1. Sensitivity analyses.',
 '(A) Adjusted CUSUM after each leave-one-out refit, plotted using original case numbers. The omitted case is shown as a gap. (B) Distribution of peak positions across 48 refits, using original case numbers. (C) Additional adjustment for ordinal POP-Q stage, coded 2, 3, or 4 as a single numeric covariate. (D) Complete centered 5-case moving averages are assigned to their center cases (cases 3–46), centered at the mean of each resulting 44-point smoothed sequence, and cumulatively summed. Cases 1–2 and 47–48 have no complete window and are not plotted; no artificial zero endpoint is added. These smoothed curves use different centering from the primary 48-case curves.')]

notes=[
 '分析对象：上传Excel的Sheet1，48例术者1机器人辅助手术。按原工作表顺序分析；日期已升序。原始上传文件未修改。用户已确认采用“总手术时间_min”列。',
 f'主要结果复现：传统峰值第31例（497.88 min）；调整后峰值第21例（488.37 min）；加入POP-Q后峰值第{cp.argmax()+1}例（{cp.max():.2f} min）。峰值是累计偏差的描述性位置，不能直接视为熟练度达标点。',
 f'10,000次置换：传统P={pr:.4f}，调整后P={pa:.4f}。正文P=0.028/0.018的精确复现需要原随机种子；本次固定种子并使用加一校正，数值差异可能来自Monte Carlo随机误差，不能为了匹配原文挑选种子。灰带为逐点95%置换参考区间，不是同时置信带。',
 f'逐例留一峰值（原病例编号）：{counts.to_dict()}。原编号与剔除后重编号不同，结果表分别列出，避免“附近”定义不清。',
 f'病例构成模型：R²={m["r2"]:.3f}，调整后R²={m["adj_r2"]:.3f}；留一交叉验证MAE={m["loo_mae"]:.2f} min。',
 'POP-Q按Aa、Ba、C、Ap、Bp中最远端的测量值及TVL计算：领先点<-1为I期，-1至+1为II期，>+1且<TVL-2为III期，≥TVL-2为IV期。本队列无I期。所得II/III/IV期为9/30/9，序号相关ρ=0.356，P=0.013，与正文一致。未使用诊断原文中的“度”直接替代POP-Q分期。',
 'POP-Q结构化测量有两处缺失：第18例Ba缺失，原文为B+3.0（点位标签需核查）；已记录Aa=+3，对本次III期判定无影响。第36例Aa缺失，原文为a+3.0，已有Ba=C=+5，TVL=7，判为IV期。未更改原表或伪造完整测量。正式发表前建议核查原病历。',
 '8组日期重复，共16例。缺少同日具体手术先后时间，当前保留工作表顺序。学习曲线依赖排序，正式发表前应核实同日先后；未对同日病例随意重排。',
 '分析时间与“手术时长_原值”有25例不一致（包含1例原值缺失），已逐例列于Data_audit；用户确认使用分析时间列，未补写或替换原值。',
 '5例居中移动平均：完整窗口+重新中心化可复现传统第32例/调整后第19例；边界缩短窗口+重新中心化为第32例/第20例。正文方法应明确选择完整窗口并以平滑后均值中心化，不应只写“5例居中移动平均”。',
 '回归采用HC3稳健协方差；95%CI及主P值采用残差自由度t分布，同时导出正态近似P值。正文应明确推断约定，不能混用不同约定的置信区间和P值。',
 '本次实际使用Python（具体版本见summary.json），没有使用R。若直接采用本次输出，论文统计软件描述应据实修改；随附脚本可复算。',
 '未生成病例筛选流程图：主表与历史表存在重复/历史信息，仅凭现有文件不能可靠推断初筛、排除人数。图1为根据研究目的新增的手术时间概览；正文明确引用的图2和补充图S1均已生成。',
 '置换局限：传统时间置换用于顺序效应的随机参考；残差置换按论文方案执行，不重新拟合模型。回归残差受杠杆值和潜在异方差影响，严格交换性未获保证；应作为探索性分析解释。短期并发症空白不能自动确认为“无”，本任务未据空白重算并发症率。',
 '交付格式：每幅图PNG 400 dpi、TIFF 600 dpi（LZW压缩）、矢量PDF和SVG；图中英文用于英文题名稿件，中文说明及可直接使用的英文图注附于报告。输出表仅包含分析序号和必要数值，不包含姓名、手术日期或病历文本。'
]
report=Document();report.add_heading('RASC学习曲线：图片与复算核查报告',0)
for t in notes:report.add_paragraph(t)
report.add_heading('手术时间趋势回归（每增加1例）',1)
table=report.add_table(rows=1,cols=4);table.style='Table Grid'
for c,t in zip(table.rows[0].cells,['模型','β (min/例)','HC3 95% CI (t)','P (t)']):c.text=t
for _,r in trends.iterrows():
    for c,t in zip(table.add_row().cells,[r['model'],f'{r.beta_per_case:.3f}',f'{r.ci_low_t:.3f} to {r.ci_high_t:.3f}',f'{r.p_t:.4f}']):c.text=t
for name,(title,caption) in zip(['Figure_1_operative_time','Figure_2_CUSUM_permutation','Figure_S1_sensitivity'],captions):
    report.add_page_break();report.add_heading(title,1);report.add_picture(str(OUT/(name+'.png')),width=Inches(6.3));report.add_paragraph(caption)
report.save(OUT/'图片与分析核查报告.docx')
(OUT/'README_图注与方法.txt').write_text('\n\n'.join(notes+['\n'+a+'\n'+b for a,b in captions])+ '\n\n复算：python analysis.py /path/to/input.xlsx /path/to/output\n依赖：numpy pandas scipy matplotlib openpyxl python-docx pillow\n',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))

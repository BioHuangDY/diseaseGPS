/**********************************************************
 *Copyright(C).2016,SJTU
 *文件：hposcan.c
 *作者：杨盛
 *版本：version 1.2.0
 *完成日期：2016/10/19 23：46
 *描述：用于辅助诊断治疗
 *更新：(1)命令行参数化
        (2)部分输出英语化
        (3)省略前向图
        (4)部分数组初始化
************************************************************/

#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include<stdlib.h>

#define Node_Num 20000 //最大HPO 所有节点数目
#define Dis_Num 9000 //最大OMIM 所有疾病数目
#define MAXLINE 10000 //每行最大长度
#define OMIMLENTH 50 //OMIM ID最大长度
#define HPLENTH 10 //HP ID最大长度
#define FILELENTH 100 //文件名最大长度
#define OPTLENTH 1000 //参数长度

typedef struct lnode{ //链表节点类型
  int No; //每个节点的Index
  double weight; //边的权值
  struct lnode* next;
} list_node, *plist_node;

typedef struct HPNode{ //HPO Term类型
  int No;
  char *id;
  char *name;
  char *alt_id;
  char *def;
  char *comment;
  char *synonym;
  char *xref;
  char *is_a;
} HPNode;

typedef struct graph_type{ //图类型
  plist_node edge[Node_Num];
  plist_node back_edge[Node_Num];
} graph_type;

HPNode hpo_term[Node_Num]; //记录所有Term
graph_type graph; //HPO 图
int visit[Dis_Num][Node_Num]; // 子图节点信息
int edges_num[Node_Num];  //HPO 每个节点的子节点数目
double eic_weight[Node_Num]; // 边的eic含量
 

typedef struct DisType{ //疾病类型
  char id[OMIMLENTH];
  char *hpidlist;
} DisType;

DisType OMIMSet[Dis_Num]; //记录所有疾病

typedef struct ScoreType{
  char id[OMIMLENTH];
  double score;
}ScoreType;


/************* 添加前向边 ****************/
/*
void AddEdge(graph_type *g,int pre_No,int pos_No,double w)
{
  plist_node temp_node=(plist_node)malloc(sizeof(list_node));
  temp_node->No=pos_No;
  temp_node->weight=w;
  temp_node->next=g->edge[pre_No];
  g->edge[pre_No]=temp_node;
}
*/
/************* 添加后向边 ****************/
void AddBackEdge(graph_type *g,int pre_No,int pos_No,double w)
{
  plist_node temp_node=(plist_node)malloc(sizeof(list_node));
  temp_node->No=pos_No;
  temp_node->weight=w;
  temp_node->next=g->back_edge[pre_No];
  g->back_edge[pre_No]=temp_node;
}

/***************** 搜索某个HP ID的index***************/
int BiSearch(int begin,int end,const char* string){
  if(strcmp(hpo_term[begin].id,string)==0)
    return begin;
  if(strcmp(hpo_term[end].id,string)==0)
    return end;
  int left=begin,right=end,mid,tmp;
  while(left+1<right){
    mid=(left+right)/2;
    tmp=strcmp(hpo_term[mid].id, string);
    if(tmp==0)
      return mid;
    else if(tmp<0)
      left=mid;
    else
      right=mid;
  }
  printf("Error: Unknown ID %s\n",string);
  exit(1);
}

/*************** 从某节点开始遍历HPO 图***************/
void DFS(const graph_type *g,int Node_No,int v[]){
  v[Node_No]=1;
  plist_node temp_node=g->back_edge[Node_No];
  while(temp_node!=NULL){
    if(v[temp_node->No]!=1)
      DFS(g,temp_node->No,v);
    temp_node=temp_node->next;
  }
}

/****************** 读取HPO文件 ********************/
int CreateHPONet(char *hpofile){
  int i;
  FILE *HPO=fopen(hpofile, "r");
  if(HPO==NULL){  //打开HPO文件失败
    printf("Error: Cannot open %s file(hpo) \n",hpofile);
    exit(1);
  }  
  char *templine=(char *)malloc(MAXLINE); //存取文件的每一行
  char *tempstr=(char *)malloc(MAXLINE); //存取行中的目标字符串
  int now=-1; //HPO term指示器
  while(fgets(templine,MAXLINE,HPO)!= NULL){
    if(sscanf(templine,"id: %s",tempstr)){ //存取id
      hpo_term[now].id=(char *)malloc(strlen(tempstr)+5);
      strcpy(hpo_term[now].id,tempstr);
    }
    if(sscanf(templine,"name: %[^\n]",tempstr)){ //存取名称
      hpo_term[now].name=(char*)malloc(strlen(tempstr)+5);
      strcpy(hpo_term[now].name,tempstr);
    }
    if(sscanf(templine,"def: %[^\n]",tempstr)){ //存取描述
      hpo_term[now].def=(char *)malloc(strlen(tempstr)+5);
      strcpy(hpo_term[now].def,tempstr);
    }
    if(sscanf(templine,"comment: %[^\n]",tempstr)){  //存取评论
      hpo_term[now].comment=(char *)malloc(strlen(tempstr)+5);
      strcpy(hpo_term[now].comment,tempstr);
    }
    if(sscanf(templine,"is_a: %s",tempstr)){ //存取父节点
      if(hpo_term[now].is_a==NULL){ //第一个父节点
	hpo_term[now].is_a=(char *)malloc(strlen(tempstr)+5);
	strcpy(hpo_term[now].is_a, tempstr);
      }       
      else{   //非首父节点，并用|隔开
	hpo_term[now].is_a=(char *)realloc(hpo_term[now].is_a, strlen(hpo_term[now].is_a)+strlen(tempstr)+5);
	strcat(hpo_term[now].is_a,"|");
	strcat(hpo_term[now].is_a,tempstr);
      }
    }

    if(strcmp(templine, "[Term]\n")==0){ //找到HPO Term
      now++;
      hpo_term[now].No=now;
      hpo_term[now].is_a=NULL;
    }
  }
  free(templine); //释放指针内存
  free(tempstr);

  int hpo_term_sum=now+1; //HPO Term总数

  /******** HPO 网络图构建 **********/
  for(i=0;i<Node_Num;i++){ //初始化网络图
    graph.edge[i]=NULL;
    graph.back_edge[i]=NULL;
  }
  int j,k=0;
  char temphp[HPLENTH]; //暂时存取每个HP ID
  int pre_no; //边的前节点
  char tempch; //暂时存取每个字符
  for(i=0;i<hpo_term_sum;i++){ //遍历所有HPO节点
    j=0;
    if(hpo_term[i].is_a!=NULL){ //忽略根节点
      while(j<=strlen(hpo_term[i].is_a)){ //按'|'截取HP ID
     	if(hpo_term[i].is_a[j]!='|' && hpo_term[i].is_a[j]!='\0'){
	  temphp[k++]=hpo_term[i].is_a[j++];
	}
	else{
	  temphp[k]='\0'; //截取成功
	  pre_no=BiSearch(0,hpo_term_sum-1,temphp); //二分查找HP ID对应的Index
	  //AddEdge(&graph,pre_no,i,1.0); //添加前向边
	  AddBackEdge(&graph,i,pre_no,1.0); //添加后向边
	  k=0;
	  j++;
	}
      }
    } 
  }
  return hpo_term_sum;
}

/***************** 读取OMIM 文件******************/
int ReadOMIM(char* omimfile){
  int i;
  FILE *OMIM=fopen(omimfile, "r");
  if(OMIM==NULL){  //打开OMIM文件
    printf("Error: Cannot Open %s file(omim)\n",omimfile);
    exit(1);
  }
  char *templine=(char *)malloc(MAXLINE); //存取文件的每一行
  char *newtempomimid=(char *)malloc(OMIMLENTH);  //存取行中的omim id 
  char *oldtempomimid=(char *)malloc(OMIMLENTH);  //存取上一行中的omim id    
  char *temphpid=(char *)malloc(HPLENTH); //存取行中的hp id 
  int now=-1; //OMIM 疾病指示器
  fgets(templine,MAXLINE,OMIM); //跳过第一行
  while(fgets(templine,MAXLINE,OMIM)!= NULL){
    sscanf(templine,"%s\t%*s\t%*s\t%s\t%*s",newtempomimid,temphpid); //截取OMIM ID 以及对应的HP ID
    if(strcmp(oldtempomimid,newtempomimid)!=0){ //读取到新的OMIM ID
      now++;
      strcpy(OMIMSet[now].id,newtempomimid);
      OMIMSet[now].hpidlist=(char *)malloc(strlen(temphpid)+5);
      strcpy(OMIMSet[now].hpidlist,temphpid);
    }
    else{
      if(strstr(OMIMSet[now].hpidlist,temphpid)==NULL){ //相同的OMIM ID,将新的HP连接到hpidlist
        OMIMSet[now].hpidlist=(char *)realloc(OMIMSet[now].hpidlist, strlen(OMIMSet[now].hpidlist)+strlen(temphpid)+5);
        strcat(OMIMSet[now].hpidlist,"|");
        strcat(OMIMSet[now].hpidlist,temphpid);
      }
    }
    strcpy(oldtempomimid,newtempomimid);
  }
  int dis_sum=now+1;
  return dis_sum;
}

/********* 初始化子图 ***************/
void SubgraphInit(int v[],int hpo_term_sum)
{
  int i;
  for(i=0;i<hpo_term_sum;i++)
    v[i]=0;
}

/************** 初始化打分数组 *************/
void ScoreInit(double v[],int dis_sum)
{
  int i;
  for(i=0;i<dis_sum;i++)
    v[i]=0;
}

/************ 所有疾病子图 ****************/
void AllSubGraph(const graph_type *g,const int hpo_term_sum,const int dis_sum )
{
  memset(visit,0,sizeof(visit));
  int i,j,k=0;
  char temphp[HPLENTH]; //存取行中的hp
  char tempch; //存取字符
  for(i=0;i<dis_sum;i++){ //遍历所有疾病
    j=0;
    while(j<=strlen(OMIMSet[i].hpidlist)){ //按'|'截取hpidlist 
      tempch=OMIMSet[i].hpidlist[j];
      if(tempch!='|' && tempch!='\0'){
	temphp[k++]=tempch;
	j++;
      }
      else{
	temphp[k]='\0'; //截取成功
      	DFS(g,BiSearch(0,hpo_term_sum-1,temphp),visit[i]); //取子图，即判断那些节点已经遍历
	k=0;
	j++;
      }
    }  
  }
}

/************* 子图的交集 *****************/
void Intersection(const int v1[],const int v2[],int v[],int hpo_term_sum)
{
  int i;
  for(i=0;i<hpo_term_sum;i++)
    if(v1[i] && v2[i])
      v[i]=1;    
}

/************* 子图的并集 *****************/
void Union(const int v1[],const int v2[],int v[],int hpo_term_sum)
{
  int i;
  for(i=0;i<hpo_term_sum;i++)
    v[i]=v1[i]+v2[i];
}

/************* HPO 每个节点的子节点数目 ************/
void CountEdges(const graph_type *g,int hpo_term_sum)
{
  int i;
  plist_node temp_node=(plist_node)malloc(sizeof(list_node));
  int count;
  for(i=0;i<hpo_term_sum;i++){
    count=0;
    if(g->back_edge[i]!=NULL){
      temp_node=g->back_edge[i];
      count++;
      while(temp_node->next!=NULL){
	temp_node=temp_node->next;
	count++;
      }
    } 
    edges_num[i]=count;
  } 
}

/***************** 子图的总边数 ********************/
int SumEdges(const int v[],int hpo_term_sum)
{
  int i;
  int sum=0;
  for(i=0;i<hpo_term_sum;i++)
    sum+=v[i]*edges_num[i]; //表型集合中所有节点的子节点数目之和
  return sum;
}

/**************** 子图间的拓扑相似性 *****************/
double TopSimilarity(const graph_type *g,const int v1[],const int v2[],int hpo_term_sum)
{
  int *v_inter=(int *)malloc(sizeof(int)*hpo_term_sum);
  SubgraphInit(v_inter,hpo_term_sum);   //子图初始化
  Intersection(v1,v2,v_inter,hpo_term_sum);    //计算v1和v2的交集，并存储在v_inter
  int s1=SumEdges(v1,hpo_term_sum);   //v1子图的总边数
  int s_inter=SumEdges(v_inter,hpo_term_sum);   //v_inter子图的总边数
  double sim_top=(double)s_inter/(double)s1;  //交集的子图总边数/v1子图总边数作为拓扑相似性
  free(v_inter);
  return sim_top;
}

/************** 边的EIC含量 ****************/
void EIC(int hpo_term_sum,int dis_sum)
{
  int i,j;
  for(j=0;j<hpo_term_sum;j++)
    for(i=0;i<dis_sum;i++)
      eic_weight[j]+=visit[i][j];
  for(i=0;i<hpo_term_sum;i++)
    if(eic_weight[i]>0)
      eic_weight[i]=-log((double)eic_weight[i]/(double)dis_sum);    //EIC为表型所包含的信息量，所注释的疾病数量除以总疾病数量的负对数
}

/***************** 子图的NIC含量 *******************/
double NIC(const graph_type *g,const int v1[],const int v2[],int hpo_term_sum)
{
  int i;
  int *v_inter=(int *)malloc(hpo_term_sum*sizeof(int));
  SubgraphInit(v_inter,hpo_term_sum);   //初始化子图
  Intersection(v1,v2,v_inter,hpo_term_sum);   //计算v1和v2的交集，并存储在v_inter   
  double NIC=0;
  for(i=0;i<hpo_term_sum;i++)
    if(v_inter[i])
      NIC+=((double)edges_num[i]*eic_weight[i]);  //表型节点的子节点数量*表型包含的信息量
  free(v_inter);
  return NIC; //NIC是所有交集表型的NIC之和
}

/****************** NIC 相似性****************/
void NICSimilarity(const graph_type *g,const int v[],double score[],int hpo_term_sum,int dis_sum)
{
  int i;
  double min_nic=10; //自己拟定
  double max_nic=0; //自己拟定
  double *nic=(double *)malloc(dis_sum*sizeof(double));
  ScoreInit(nic,dis_sum); //初始化打分数组
  for(i=0;i<dis_sum;i++){
    nic[i]=NIC(g,v,visit[i],hpo_term_sum);  //设定nic的上下限。NIC相似性是疾病与表型之间的相似性，RealtiveBestPair
    max_nic=(nic[i]>max_nic)?nic[i]:max_nic;
    min_nic=(nic[i]<min_nic)?nic[i]:min_nic;
  }
  for(i=0;i<dis_sum;i++)
    score[i]=(nic[i]-min_nic)/(max_nic-min_nic);  //归一化，每一种疾病自行进行归一化
  free(nic);
}

/************** 子图的MEIC含量 *****************/
double MEIC(const graph_type *g,const int v1[],const int v2[],int hpo_term_sum)
{
  int i;
  int *v_inter=(int *)malloc(hpo_term_sum*sizeof(int));
  SubgraphInit(v_inter,hpo_term_sum);
  double MEIC=0;
  double tempeic=0;
  Intersection(v1,v2,v_inter,hpo_term_sum);
  for(i=0;i<hpo_term_sum;i++)
    if(v_inter[i]){
      tempeic=eic_weight[i];
      MEIC=(tempeic>MEIC)?tempeic:MEIC; //两个表型集合的交集部分，信息量最多的作为两者的MEIC
    }
  free(v_inter);
  return MEIC;
}

/****************** MEIC 相似性****************/
void MEICSimilarity(const graph_type *g,const int v[],double score[],int hpo_term_sum,int dis_sum)
{
  int i;
  double min_meic=10; //自己拟定
  double max_meic=0; //自己拟定
  double *meic=(double *)malloc(dis_sum*sizeof(double));
  ScoreInit(meic,dis_sum);
  for(i=0;i<dis_sum;i++){
    meic[i]=MEIC(g,v,visit[i],hpo_term_sum);
    max_meic=(meic[i]>max_meic)?meic[i]:max_meic;
    min_meic=(meic[i]<min_meic)?meic[i]:min_meic;
  }
  for(i=0;i<dis_sum;i++)
    score[i]=(meic[i]-min_meic)/(max_meic-min_meic);  //表型集合与疾病之间的MEIC相似性使用表型注释中信息量最高的来表示
  free(meic);
}

/****************  查询 ****************/
void Query(const graph_type *g,int v[],int hpo_term_sum ,char *hpidlist)
{
  int i=0,j=0;
  char tempch;
  char *temphp=(char *)malloc(HPLENTH*sizeof(char));
  while(i<=strlen(hpidlist)){ //按','截取hpidlist
    tempch=hpidlist[i];
    if(tempch!=',' && tempch!='\0'){
	temphp[j++]=tempch;
	i++;
      }
      else{
	temphp[j]='\0'; //截取成功
	DFS(g,BiSearch(0,hpo_term_sum-1,temphp),v);
	j=0;
	i++;
      }
  }
}


/*************** 排序比较 *************/
int Compare(const void*p1,const void *p2)
{
  return(*(ScoreType*)p2).score>(*(ScoreType*)p1).score?1:-1;
}

int main(int argc,char *argv[]){

  /*********** 处理命令行参数 ****************/
  int i,optflag=0,helpflag=0;
  char *tempopt=(char *)malloc(OPTLENTH*sizeof(char));
  char *omimfile=(char *)malloc(OPTLENTH*sizeof(char));
  char *hpofile=(char *)malloc(OPTLENTH*sizeof(char));
  char *hpidlist=(char *)malloc(OPTLENTH*sizeof(char));
  char *outfile=(char *)malloc(OPTLENTH*sizeof(char));
  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-help")==0){
      helpflag=1;
      printf("Usage    : ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
      printf(" -omim   : Input OMIM file name\n");
      printf(" -obo    : Input HPO file name\n");
      printf(" -term   : Input HPO terms you want to query\n");
      printf(" examples: HP:0000001,HP:000002,HP:0000003\n");
      printf(" -o      : Input output file\n");
      printf(" -h/-help: Usage for program\n");
      return 1;
    }
  }
  if(argc!=7 && argc!=9 && helpflag==0){
    printf("Error:  Too many or too lees arguements\n");
    printf("Usage:  ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
    return 1;
  }
  for(i=1;i<argc;i++){
    if(strstr(argv[i],"-")){
      sscanf(argv[i],"-%s",tempopt);
      optflag=1;
      if(strcmp(tempopt,"omim")==0){
	if((i!=argc-1) && argv[i+1][0]!='-'){
	  strcpy(omimfile,argv[i+1]);
	  i+=1;
	}
	else{
	  printf("Error: No input OMIM file\n");
	  printf("Usage:  ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
	  return 1;
	}
      }
      else if(strcmp(tempopt,"obo")==0){
	if((i!=argc-1) && argv[i+1][0]!='-'){
	  strcpy(hpofile,argv[i+1]);
	  i+=1;
	}
	else{
	  printf("Error: No input HPO file\n");
	  printf("Usage:  ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
	  return 1;
	}
      }
      else if(strcmp(tempopt,"term")==0){
	if((i!=argc-1) && argv[i+1][0]!='-'){
	  strcpy(hpidlist,argv[i+1]);
	  i+=1;
	}
	else{
	  printf("Error: No input query terms\n");
	  printf("Usage:  ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
	  return 1;
	}
      }
      else if(strcmp(tempopt,"o")==0){
	if((i!=argc-1) && argv[i+1][0]!='-'){
	  strcpy(outfile,argv[i+1]);
	  i+=1;
	}
	else{
	  printf("Error: No input outputfile\n");
	  printf("Usage:  ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
	  return 1;
	}
      }
      else if(strcmp(tempopt,"help")==0 || strcmp(tempopt,"h")==0){
	printf("Usage    : ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
	printf(" -omim   : Input OMIM file name\n");
	printf(" -obo    : Input HPO file name\n");
	printf(" -term   : Input HPO terms you want to query\n");
	printf(" examples: HP:0000001,HP:000002,HP:0000003\n");
	printf(" -o      : Input output file\n");
	printf(" -h/-help: Usage for program\n");
	return 1;
      }
      else{
	printf("Error: Unknown options -%s\n",tempopt);
	printf("Usage:  ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
	return 1;
      }
    }
  }
  if(!optflag){
    printf("Error: No options\n");
    printf("Usage:  ./hpo -omim omimfile -obo hpofile -term hpterms -o outfile\n");
    return 1;
  }
  

  /***************** 建立网络以及分析 ******************/
  int hpo_term_sum=CreateHPONet(hpofile);
  int dis_sum=ReadOMIM(omimfile);
  //printf("%d\n",dis_sum);
  AllSubGraph(&graph,hpo_term_sum,dis_sum);
  CountEdges(&graph,hpo_term_sum);
  EIC(hpo_term_sum,dis_sum);
  double b=0;
  double *meic_score=(double *)malloc(dis_sum*sizeof(double));
  double *top_score=(double *)malloc(dis_sum*sizeof(double));
  ScoreType *balance_score=(ScoreType *)malloc(dis_sum*sizeof(ScoreType));
  int *query_v=(int *)malloc(hpo_term_sum*sizeof(int));
  SubgraphInit(query_v,hpo_term_sum);
  ScoreInit(meic_score,dis_sum);
  ScoreInit(top_score,dis_sum);
  double coeff=0.9;
  Query(&graph,query_v,hpo_term_sum,hpidlist);
  MEICSimilarity(&graph,query_v,meic_score,hpo_term_sum,dis_sum);
  for(i=0;i<dis_sum;i++){
    top_score[i]=TopSimilarity(&graph,query_v,visit[i],hpo_term_sum);
    strcpy(balance_score[i].id,OMIMSet[i].id);
    balance_score[i].score=coeff*top_score[i]+(1-coeff)*meic_score[i];  //拓扑相似性与MEIC相似性以9：1的比例相加
  }
  qsort(balance_score,dis_sum,sizeof(balance_score[0]),Compare);
  if(strcmp(outfile,"")!=0){
    FILE *OUT=fopen(outfile,"w");
    for(i=0;i<dis_sum;i++){
      fprintf(OUT,"%s\t%lf\n",balance_score[i].id,balance_score[i].score);
    }
    fclose(OUT);
  }
  else{
    for(i=0;i<dis_sum;i++)
      printf("%s\t%lf\n",balance_score[i].id,balance_score[i].score);
  }
  //printf("**************** Analysis Completed *********************\n");
  return 0;
}

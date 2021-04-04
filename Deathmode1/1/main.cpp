#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <GLUT/glut.h>
#include <time.h>
#include "Monitor.hpp"
//#include <opencv/cv.h>
#include "Vector.h"
#include <string.h>
//#include "DefineCell.h"
#include <unistd.h>
#include <sys/stat.h>
//#include <windows.h>

#include <png.h>

int WinID[3] = {};
double dt = 0.0001;      //dt = 0.01;
#define MAX_STEP (100000) //500000最大タイムステップ数700000
#define N (2)            //一辺のマス数,N*N＝セル数
//#define N0 (100)
void capture(int *);
//ディスプレイ
int displayInterval = 50000;

#define zoom2 (66) //数値大→遠くなる
//#define xcenter1 (100.0)         //大→左へずれる（フィールドが）
//#define ycenter1 (80.0)          //大→下へずれる
#define xcenter2 (170.0) //大→左へずれる（フィールドが）
#define ycenter2 (100.0) //大→下へずれる

Monitor monitor;

//要素
double u_old[N] = {};
double v_old[N] = {};
double u_new[N] = {};
double v_new[N] = {};

double eps = 0.001, au = 37.0, c = 5.0, sigma0 = 0.0, av = 0.1 * au, mu0 = 0.01; //mu0: Coupling Strength
double u0 = 10.0, v0 = 0.001;
double vCoupling = 0.3;

int csv_read = 0;
int FlagDataSave = 1;

int tsJump = 0;
int Flaaaag = 0, count1 = 0, count2 = 0;

int ParameterSetCount = 0, EndSetCount = 48;

double t;
int ts = 0, pictureCount = 0; //現在のタイムステップ

void mouse(int, int, int, int);
void init();
void idle();

FILE *fpValues, *fpV, *fpU, *fp, *fpZ, *fpCellcount, *fpIbhibitedCellcount;
char filename[256] = "", strFolderParaset[150] = "", strParaSetTxt[256] = "", strCellCount[256] = "";
int FileCount = 0;

void ParameterInitialization()
{
  char tmp[256];
  char *p, *p2, *r;
  char **endptr = NULL;
  int fileNum = 0, i = 0;
  char strPos[128];

  //現在のディレクトリを取得
  p = getcwd(tmp, sizeof(tmp) - 1);
  printf("tmp=%s\n", tmp);
  strPos[0] = '/';
  strcpy(strPos + 1, tmp);

  //フォルダの番号を探す
  char str[] = "";
  int sizeNum = 0, sizeCharas = 0;
  r = strtok(strPos, "/"); //一文字目を切り出します。
  sizeCharas = strlen(p);
  strcat(str, "/");
  strcat(str, r);
  while (r = strtok(NULL, "/"))
  { //次の文字を切り出します。無ければNULLを返します。
    fileNum = (int)(strtol(r, endptr, 10));
    sizeNum = strlen(r);
    printf("%s\n", r);
    strcat(str, "/");
    strcat(str, r);
  }
  printf("fileNum=%d\n", fileNum);
  printf("%s\n", str);
  printf("sizeCharas=%d\n", sizeCharas);
  printf("sizeNum=%d\n", sizeNum);
  printf("Check00\n");
  //parameters.csvのディレクトリを抜き出す
  char strDir[sizeCharas - sizeNum + 1];
  strncpy(strDir, tmp, sizeCharas - sizeNum);
  sizeCharas - sizeNum + 1;
  strDir[sizeCharas - sizeNum] = '\0';
  printf("strDir=%s\n", strDir);
  printf("Check3\n");
  //parameterset.csvをオープンする
  strcat(strDir, "/parameters.csv");
  char *fname = strDir; //"/Users/mikamihiroshi/Document/Research/programs/check/parameters.csv";
  FILE *fp;
  /*
  fp = fopen(fname, "r");
  if (fp == NULL)
  {
    printf("%sファイルが開けません\n", fname);
    exit(0);
  }
  printf("Check1\n");
  //csvファイルからパラメータ値を代入する
  char temp;
  char StrCs1[200], StrCs[200], StrCv[200], StrCs2[200];
  double kDz = 0.0;
  fscanf(fp, "%[^,],%s", StrCs, StrCs1, StrCs2, StrCv);
  printf("Check2\n");

  for (i = 0; i < fileNum; i++)
  { //csvファイルの1行目の数値はなぜか読み取ってくれない
    fscanf(fp, "%lf,%d", &kDz, &NInitRadius);
  }*/
  //cs = strtod(StrCs, endptr);

  if (Flaaaag == 0)
  { //シミュレーション開始時の処理
    Flaaaag++;
    ///パラメータセットごとに数値を変更するパラメータ(2つ)
    count1 = ParameterSetCount % 6;
    count2 = ParameterSetCount / 6;
    //printf("0 count1=%d count2=%d\n", count1, count2);
  }
  if (Flaaaag >= 2)
  { //パラメータの更新
    FileCount = 0;
    //countの更新
    if (count1 % 5 == 0 && count1 != 0)
    {
      count1 = 0;
      count2++;
    }
    else
    {
      count1++;
    }
    //printf("ParameterSetCount=%d", ParameterSetCount);
    //printf("count1=%d count2=%d\n", count1, count2);
  }
  else
  {
    Flaaaag++;
  }

  //パラメータの更新
  double vCouplingArray[8] = {0.0001,0.001,0.01,0.1,1.0,10.0,100.0,1000.0};
  double mu0Array[6] = {0.0001,0.001,0.01,0.1,1.0,10.0};
  vCoupling = vCouplingArray[count2];
  mu0 = mu0Array[count1];
  printf("Set=%d  count1=%d count2=%d vCoupling=%lf  mu0=%lf\n",ParameterSetCount,count1,count2,vCoupling,mu0);
}

void FolderAndFileInit()
{
  char tmp[256] = "";
  char *p;
  int i = 0, j = 0;
  p = getcwd(tmp, sizeof(tmp) - 1);

  if (p == NULL)
  {
    printf("ERROR,getcwd() ret=NULL\n");
    perror("getcwd\n");
    exit(EXIT_FAILURE);
  }

  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  if (mkdir(strFolderParaset, S_IRWXU) == 0)
  { //フォルダ作成
    printf("フォルダ作成に成功しました。\n");
  }
  else
  {
    do
    {
      ParameterSetCount++;
      sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
    } while (mkdir(strFolderParaset, S_IRWXU) != 0);
    /* printf("フォルダ作成に失敗しました。\n");
    printf("%s\n", strFolderParaset);
    exit(0);*/
  }
  ParameterInitialization();

  //パラメータ値を保存
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(strParaSetTxt, "%s/ParaSet.txt", strFolderParaset);
  fpValues = fopen(strParaSetTxt, "w");
  if (fpValues == NULL)
  {
    printf("ファイル作成失敗fpV\n");
    printf("filename=%s\n", strParaSetTxt);
    exit(0);
  }
  fprintf(fpValues, "eps=%lf,au=%lf,c=%lf,sigma0=%lf,av=%lf,mu0=%lf\n", eps, au, c, sigma0, av, mu0);
  fclose(fpValues);
  /////
  sprintf(strFolderParaset, "%s/Results/ParameterSet%d", tmp, ParameterSetCount);
  sprintf(strCellCount, "%s/Cellcounts%d.csv", strFolderParaset, ParameterSetCount);
  fpCellcount = fopen(strCellCount, "w");
  if (fpCellcount == NULL)
  {
    printf("ファイル作成失敗fpV\n");
    printf("filename=%s\n", strCellCount);
    exit(0);
  }
  if (N == 2)
  {
    fprintf(fpCellcount, "%d,%lf,%lf,%lf,%lf\n", ts, u_new[0], u_new[1], v_new[0], v_new[1]);
  }
  else
  {
    fprintf(fpCellcount, "%d,%lf,%lf\n", ts, u_new[0], v_new[0]);
  }
  fclose(fpCellcount);
  ////////////////
}

void init()
{
  //変数の無次元化
  int x = 0;
  ts = 0;
  t = 0.0;
  srand(time(NULL));
  

  if (displayInterval > MAX_STEP)
  {
    displayInterval = MAX_STEP;
  }

  for (x = 0; x < N; x++)
  {
    /*//vだけ振動した au=av=37で
    u_old[x] = 0.000010, u_new[x] = u_old[x];
    v_old[x] = 1.0, v_new[x] = v_old[x];*/
    u0 = 9.0 + (rand() * (11.0 - 9.0 + 1.0) / (1.0 + RAND_MAX));
    v0 = (5.0 + (rand() * (15.0 - 5.0 + 1.0) / (1.0 + RAND_MAX))) * 0.001;
    u_old[x] = u0, u_new[x] = u_old[x];
    u0 = 9.0 + (rand() * (11.0 - 9.0 + 1.0) / (1.0 + RAND_MAX));
    v0 = (5.0 + (rand() * (15.0 - 5.0 + 1.0) / (1.0 + RAND_MAX))) * 0.001;
    v_old[x] = v0, v_new[x] = v_old[x];
  }

  if (FlagDataSave == 1)
  {
    FolderAndFileInit();
  }
}

void FieldCycle()
{

  //拡散項を求める
  int x = 0;
  double Expornential = 0.0;
  double FirstAndSecondTerm_u = 0.0, ThirdTerm_u = 0.0, ForthTerm_u = 0.0;
  //double ThirdTerm_v = 0.0;
  for (x = 0; x < N; x++)
  {
    Expornential = exp(u_old[x] / (1 + (u_old[x] / c)));
    FirstAndSecondTerm_u = (1.0 / eps) * (-u_old[x] + au * v_old[x] * Expornential);
    ThirdTerm_u = sigma0 * pow(1.0 + u_old[x] / c, 4.0);
    ForthTerm_u = sigma0 * mu0 * pow(1.0 + u_old[1 - x] / c, 4.0);
    if (N == 2)
    {
      u_new[x] += (FirstAndSecondTerm_u - ThirdTerm_u + ForthTerm_u) * dt;                                      //Heat
      v_new[x] += (1.0 - v_old[x] - av * v_old[x] * Expornential + vCoupling * (v_old[1 - x]-v_old[x])) * dt; //Air
    }
    else
    {
      u_new[x] += (FirstAndSecondTerm_u - ThirdTerm_u) * dt;
      v_new[x] += (1.0 - v_old[x] - av * v_old[x] * Expornential) * dt;
    }
  }

  for (x = 0; x < N; x++)
  {
    u_old[x] = u_new[x];
    v_old[x] = v_new[x];
  }
  ts++;
  pictureCount++;
}

void keyboard(unsigned char key, int x, int y)
{
  double tmpcin;
  switch (key)
  {
  case 'q':
    exit(0);
    break;
  case '\033': /* '\033' = ESC */
    exit(0);
    break;
  case 'h':
    monitor.SetCenter(1.0, 0);
    break;
  case 'l':
    monitor.SetCenter(-1.0, 0);
    break;
  case 'j':
    monitor.SetCenter(0, 10.0);
    break;
  case 'k':
    monitor.SetCenter(0, -10.0);
    break;
  case 'g':
    break;
  case 'z':
    monitor.SetZoom(1.1 / 1.0);
    break;
  case 'x':
    monitor.SetZoom(1.0 / 1.1);
    break;
  }
}

void idle(void)
{
  int i = 0, j = 0;


    FieldCycle(); //シミュレーションを回す
    fpCellcount = fopen(strCellCount, "a");
    if (fpCellcount == NULL)
    {
      printf("ファイル作成失敗fpV\n");
      printf("filename=%s\n", strCellCount);
      exit(0);
    }
    if (N == 2)
    {
      fprintf(fpCellcount, "%d,%lf,%lf,%lf,%lf\n", ts, u_new[0], u_new[1], v_new[0], v_new[1]);
    }
    else
    {
      fprintf(fpCellcount, "%d,%lf,%lf\n", ts, u_new[0], v_new[0]);
    }

    fclose(fpCellcount);
 

  if (ts >= MAX_STEP)
  {
    ParameterSetCount++;
    if (ParameterSetCount >= EndSetCount)
    {
      exit(0);
    }
    else
    {
      init();
    }
  }

  if (ts >= tsJump)
  {
    for (i = 0; i < 0; i++)
    {
      glutSetWindow(WinID[i]);
      glutPostRedisplay();
    }
  }
}

void display1(void)
{
}

void display2(void)
{
}

void display3(void)
{
}

void mouse(int button, int state, int x, int y)
{
  switch (button)
  {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
    }
    else
    {
      glutIdleFunc(idle);
    }
    break;

  case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
      std::cout << "middle: on" << std::endl;
    }
    else
    {
      glutIdleFunc(idle);
      std::cout << "middle: off" << std::endl;
    }
    break;

  case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN)
    {
      glutIdleFunc(0);
      std::cout << "right: on" << std::endl;
    }
    else
    {
      glutIdleFunc(idle);
      std::cout << "right: off" << std::endl;
    }
    break;
  }
}

void resize1(int w, int h)
{
  glViewport(0, 0, w, h);
  //glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-w / 0.0, w / 20.0, -h / 20.0, h / 20.0, -3.0, 3.0);
  //gluPerspective( 30.0, (double)w / (double)h, 1.0 , 100.0 );
  //gluLookAt( 0.0 , 0.0 , 3.8 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 );
  //monitor.SetWindowSize( w , h );
  //glMatrixMode(GL_MODELVIEW);
}

void resize2(int w, int h)
{
  glViewport(0, 0, w, h);

  //glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-w / 0.0, w / 10.0, -h / 20.0, h / 20.0, -3.0, 3.0);
  //gluPerspective( 30.0, (double)w / (double)h, 1.0 , 100.0 );
  //gluLookAt( 0.0 , 0.0 , 3.8 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 );
  //monitor.SetWindowSize( w , h );
  //glMatrixMode(GL_MODELVIEW);
}

void OpenGL_init(int *argcp, char **argv)
{
  init(); //初期条件を設定
          /*
  //ウィンドウ２つ目
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  monitor.SetWindowSize(400, 300);
  monitor.SetMovieMode(1);
  monitor.SetZoom(zoom2);
  monitor.SetCenter(xcenter2, ycenter2);
  glutInitWindowSize(monitor.GetWindowSize(Monitor::X), monitor.GetWindowSize(Monitor::Y));
  glutInitWindowPosition(600, 100);
  WinID[1] = glutCreateWindow("simulationX");
  glutDisplayFunc(display2);
  glutReshapeFunc(resize1);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glClearColor(1.0, 1.0, 1.0, 1.0);

  //ウィンドウ２つ目
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  monitor.SetWindowSize(400, 300);
  monitor.SetMovieMode(1);
  glutInitWindowSize(monitor.GetWindowSize(Monitor::X), monitor.GetWindowSize(Monitor::Y));
  glutInitWindowPosition(1200, 100);
  WinID[2] = glutCreateWindow("simulationY");
  glutDisplayFunc(display3);
  glutReshapeFunc(resize1);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glClearColor(1.0, 1.0, 1.0, 1.0);
*/
  //ウィンドウ１つめ（これを下にしたら保存できるようになった）
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  /*monitor.SetWindowSize(400, 300);
  monitor.SetMovieMode(1);
  monitor.SetZoom(zoom1);
  monitor.SetCenter(xcenter1, ycenter1);*/
  glutInitWindowSize(monitor.GetWindowSize(Monitor::X), monitor.GetWindowSize(Monitor::Y));
  glutInitWindowPosition(10, 100);
  WinID[0] = glutCreateWindow("simulation");
  glutDisplayFunc(display1);
  glutReshapeFunc(resize2);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glClearColor(1.0, 1.0, 1.0, 1.0);
}

/*void monitor_init()
{
  monitor.SetWindowSize(400, 300);
  monitor.SetMovieMode(1);
  monitor.SetZoom(zoom);
  monitor.SetCenter(xcenter, ycenter);
}*/

void capture(int *pts)
{
  char filepath[100]; //= "./MovieDir/output.png";
  sprintf(filepath, "./MovieDir/%d.png", *pts);
  png_bytep raw1D;
  png_bytepp raw2D;
  int i;
  int width = glutGet(GLUT_WINDOW_WIDTH);
  int height = glutGet(GLUT_WINDOW_HEIGHT);

  // 構造体確保
  FILE *fp = fopen(filepath, "wb");
  png_structp pp = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  png_infop ip = png_create_info_struct(pp);
  // 書き込み準備
  png_init_io(pp, fp);
  png_set_IHDR(pp, ip, width, height,
               8,                   // 8bit以外にするなら変える
               PNG_COLOR_TYPE_RGBA, // RGBA以外にするなら変える
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // ピクセル領域確保
  raw1D = (png_bytep)malloc(height * png_get_rowbytes(pp, ip));
  raw2D = (png_bytepp)malloc(height * sizeof(png_bytep));
  for (i = 0; i < height; i++)
    raw2D[i] = &raw1D[i * png_get_rowbytes(pp, ip)];
  // 画像のキャプチャ
  glReadBuffer(GL_FRONT);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 初期値は4
  glReadPixels(0, 0, width, height,
               GL_RGBA,          // RGBA以外にするなら変える
               GL_UNSIGNED_BYTE, // 8bit以外にするなら変える
               (void *)raw1D);
  // 上下反転
  for (i = 0; i < height / 2; i++)
  {
    png_bytep swp = raw2D[i];
    raw2D[i] = raw2D[height - i - 1];
    raw2D[height - i - 1] = swp;
  }
  // 書き込み
  png_write_info(pp, ip);
  png_write_image(pp, raw2D);
  png_write_end(pp, ip);
  // 開放
  png_destroy_write_struct(&pp, &ip);
  fclose(fp);
  free(raw1D);
  free(raw2D);

  //printf("write out screen capture to '%s'\n", filepath);
}

int main(int argc, char *argv[])
{
  //monitor_init();
  std::cout << "monitor init OK" << std::endl;
  glutInit(&argc, argv);

  OpenGL_init(&argc, argv);
  std::cout << "OpenGL init OK" << std::endl;

  // glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);
  glutMainLoop(); //無限ループ

  return 0;
}

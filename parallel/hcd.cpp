/****************************************************************************
**
** For Copyright & Licensing information, see COPYRIGHT in project root
**
****************************************************************************/

#include "hcd.h"

#include <QVector>
#include <QByteArray>
#include <QColor>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <iostream>

#include "Timer.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

typedef int v8si __attribute__ ((vector_size (32)));
static v8si one = { 1, -1, 1, -1, 1, -1, 1, -1 };

/****************************************************************************
                __   ___                 __  __           __
     ___  __ __/ /  / (_)___  __ _  ___ / /_/ /  ___  ___/ /__
    / _ \/ // / _ \/ / / __/ /  ' \/ -_) __/ _ \/ _ \/ _  (_-<
   / .__/\_,_/_.__/_/_/\__/ /_/_/_/\__/\__/_//_/\___/\_,_/___/
  /_/

****************************************************************************/


/****************************************************************************
**
** Author: Marc Bowes
**
** Detects circles in the specified QImage
**
****************************************************************************/
QImage HoughCircleDetector::detect(const QImage &source, unsigned int min_r, unsigned int max_r)
{
  QImage binary = edges(source);
  QImage detection = source.convertToFormat(QImage::Format_RGB888);
  
  /* build a vector to hold images in Hough-space for radius 1..max_r, where
  max_r is specified or the maximum radius of a circle in this image */
  if(min_r == 0) min_r = 5;
  if(max_r == 0) max_r = MIN(source.width(), source.height()) / 2;

  /* find all the edges */
  PointArray edge;
  for(unsigned int y = 0; y < binary.height(); y++)
    for(unsigned int x = 0; x < binary.width(); x++)
      if(binary.pixelIndex(x, y) == 1)
        edge.append(QPoint(x, y));

  QVector<Image> houghs(max_r - min_r);
  
  #pragma omp parallel for
  for(unsigned int i = min_r; i < max_r; i++)
  {
#if TIMER
    ggc::Timer t("radius");
    t.start();
    unsigned int tid = omp_get_thread_num();
#endif

    /* instantiate Hough-space for circles of radius i */
    Image &hough = houghs[i - min_r];
    hough.resize(binary.width() + 2 * i);
    for(unsigned int x = 0; x < hough.size(); x++)
    {
      hough[x].resize(binary.height() + 2 * i);
      hough[x].fill(0);
    }
    
    /* unroll accumulate circle */
    PointArray::const_iterator iter = edge.constBegin();
    if (edge.size() % 2 == 1)
        accum_circle(hough, *iter++, i);
    for (; iter != edge.constEnd();)
        accum_circle2(hough, *iter++, *iter++, i);
    
    /* loop through all the Hough-space images, searching for bright spots, which
    indicate the center of a circle, then draw circles in image-space */
    unsigned int threshold = 4.9 * i;
    for(unsigned int x = i; x < hough.size() - i; x++)
    {
      for(unsigned int y = i; y < hough[x].size() - i; y++)
      {
        if(hough[x][y] > threshold)
        {
          draw_circle(detection, QPoint(x - i, y - i), i, Qt::yellow);
        }
      }
    }
#if TIMER
    t.stop();
    printf("THREAD %u Radius %d Time: %llu ns\n", tid, i, t.duration());
#endif
  }
    
  return detection;
}


/****************************************************************************
               _           __                  __  __           __
     ___  ____(_)  _____ _/ /____   __ _  ___ / /_/ /  ___  ___/ /__
    / _ \/ __/ / |/ / _ `/ __/ -_) /  ' \/ -_) __/ _ \/ _ \/ _  (_-<
   / .__/_/ /_/|___/\_,_/\__/\__/ /_/_/_/\__/\__/_//_/\___/\_,_/___/
  /_/

****************************************************************************/

/****************************************************************************
**
** Author: Marc Bowes
**
** Accumulates a circle on the specified image at the specified position with
** the specified radius, using the midpoint circle drawing algorithm
**
** Adapted from: http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
**
****************************************************************************/
void HoughCircleDetector::accum_circle(Image &image, const QPoint &position, unsigned int radius)
{
  int cx = position.x() + radius;
  int cy = position.y() + radius;

  int r = radius;
  int f = 1 - r;
  int ddF_x = 1;
  int ddF_y = -2 * r;
  int x = 0;
  int y = r;

  image[cx][cy + r]++;
  image[cx][cy - r]++;
  image[cx + r][cy]++;
  image[cx - r][cy]++;
  
  while(x < y)
  {
    if(f >= 0) { y--; ddF_y += 2; f += ddF_y; }
    x++; ddF_x += 2; f += ddF_x; 

    image[cx + x][cy + y]++;
    image[cx - x][cy + y]++;
    image[cx + x][cy - y]++;
    image[cx - x][cy - y]++;
    image[cx + y][cy + x]++; 
    image[cx + y][cy - x]++;
    image[cx - y][cy + x]++; 
    image[cx - y][cy - x]++;
  }
}

void HoughCircleDetector::accum_circle2(Image &image, const QPoint &p1, const QPoint &p2, unsigned int radius)
{
  int cx1 = p1.x() + radius; int cx2 = p2.x() + radius;
  int cy1 = p1.y() + radius; int cy2 = p2.y() + radius;

  int r = radius;
  int f = 1 - r;
  int ddF_x = 1;
  int ddF_y = -2 * r;
  int x = 0;
  int y = r;

  image[cx1][cy1 + r]++; image[cx2][cy2 + r]++;
  image[cx1][cy1 - r]++; image[cx2][cy2 - r]++;
  image[cx1 + r][cy1]++; image[cx2 + r][cy2]++;
  image[cx1 - r][cy1]++; image[cx2 - r][cy2]++;
  
  while(x < y)
  {
    if(f >= 0) { y--; ddF_y += 2; f += ddF_y; }
    x++; ddF_x += 2; f += ddF_x; 

    image[cx1 + x][cy1 + y]++; image[cx2 + x][cy2 + y]++;
    image[cx1 - x][cy1 + y]++; image[cx2 - x][cy2 + y]++;
    image[cx1 + x][cy1 - y]++; image[cx2 + x][cy2 - y]++;
    image[cx1 - x][cy1 - y]++; image[cx2 - x][cy2 - y]++;
    image[cx1 + y][cy1 + x]++; image[cx2 + y][cy2 + x]++;  
    image[cx1 + y][cy1 - x]++; image[cx2 + y][cy2 - x]++;
    image[cx1 - y][cy1 + x]++; image[cx2 - y][cy2 + x]++;  
    image[cx1 - y][cy1 - x]++; image[cx2 - y][cy2 - x]++;
  }
}

/****************************************************************************
**
** Author: Marc Bowes
**
** Draws a circle on the specified image at the specified position with
** the specified radius, using the midpoint circle drawing algorithm
**
** Adapted from: http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
**
****************************************************************************/
void HoughCircleDetector::draw_circle(QImage &image, const QPoint &position, unsigned int radius, const QColor &color)
{
  QRgb rgb = color.rgb();

  int cx = position.x();
  int cy = position.y();

  int f = 1 - radius;
  int ddF_x = 1;
  int ddF_y = -2 * radius;
  int x = 0;
  int y = radius;
  
  draw_pixel(image, cx + radius, cy, rgb);
  draw_pixel(image, cx - radius, cy, rgb);
  draw_pixel(image, cx, cy + radius, rgb);
  draw_pixel(image, cx, cy - radius, rgb);
  
  while(x < y)
  {
    if(f >= 0) { y--; ddF_y += 2; f += ddF_y; }
    x++; ddF_x += 2; f += ddF_x;

    draw_pixel(image, cx + x, cy + y, rgb);
    draw_pixel(image, cx - x, cy + y, rgb);
    draw_pixel(image, cx + x, cy - y, rgb);
    draw_pixel(image, cx - x, cy - y, rgb);
    draw_pixel(image, cx + y, cy + x, rgb);
    draw_pixel(image, cx + y, cy - x, rgb);
    draw_pixel(image, cx - y, cy + x, rgb);
    draw_pixel(image, cx - y, cy - x, rgb);
  }
}

/****************************************************************************
**
** Author: Marc Bowes
**
** Draws at the specified position
**
****************************************************************************/
void HoughCircleDetector::draw_pixel(QImage &image, const int x, const int y, const QRgb &rgb)
{
  /* bounds checking */
  if(x < 0 || x >= image.width() || y < 0 || y >= image.height())
  { return; }
  
  image.setPixel(x, y, rgb);
}

/****************************************************************************
**
** Author: Marc Bowes
**
** Detects edges in the specified QImage
**
****************************************************************************/
QImage HoughCircleDetector::edges(const QImage &source)
{
  /* initialisation */
  QImage binary = QImage(source.size(), QImage::Format_Mono);
  
  /*** Sobel edge detection ***/
  
  /* Set up Lx, Ly */
  QVector<QByteArray> Lx(3), Ly(3);
  
  Lx[0][0] = -1;  Lx[0][1] = +0;  Lx[0][2] = +1;
  Lx[1][0] = -2;  Lx[1][1] = +0;  Lx[1][2] = +2;
  Lx[2][0] = -1;  Lx[2][1] = +0;  Lx[2][2] = +1;
  
  Ly[0][0] = +1;  Ly[0][1] = +2;  Ly[0][2] = +1;
  Ly[1][0] = +0;  Ly[1][1] = +0;  Ly[1][2] = +0;
  Ly[2][0] = -1;  Ly[2][1] = -2;  Ly[2][2] = -1;

#if TIMER
    ggc::Timer t("sobel");
    t.start();
    unsigned int tid = omp_get_thread_num();
#endif
  
  #pragma omp parallel for
  for(unsigned int y = 0; y < source.height(); y++)
  {
    for(unsigned int x = 0; x < source.width(); x++)
    {
      double new_x = 0.0, new_y = 0.0;
      
      /* gradient */
      for(int i = -1; i <= 1; i++)
      {
        for(int j = -1; j <= 1; j++)
        {
          /* these are offset co-ords */
          int _x = x + i;
          int _y = y + j;
          
          /* bounds checking */
          if (_x < 0)                     _x = -_x;
          else if (_x >= source.width())  _x = 2 * source.width() - _x - 2;
          
          if (_y < 0)                     _y = -_y;
          else if (_y >= source.height()) _y = 2 * source.height() - _y - 2;
          
          /* accumulate */
          int gray = qGray(source.pixel(_x, _y));
          new_x += Lx[i + 1][j + 1] * gray;
          new_y += Ly[i + 1][j + 1] * gray;
        }
      }
      
      /* using 128 as a threshold, decide if the steepness is sufficient (= edge = 1) */
      int pixel = sqrt(pow(new_x, 2) + pow(new_y, 2)) > 128 ? 1 : 0;
      binary.setPixel(x, y, pixel);
    }
  }

#if TIMER
    t.stop();
    printf("Sobel Time: %llu ns\n", t.duration());
#endif
  
  return binary;
}


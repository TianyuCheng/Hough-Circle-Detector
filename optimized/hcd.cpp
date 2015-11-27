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

// #define GCC_5_2
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

typedef int v8si __attribute__ ((vector_size (32)));

static v8si *lower_bounds;
static v8si *higher_bounds;

bool QPointLessThan(const QPoint &s1, const QPoint &s2)
{
    return s1.x() == s2.x() ? s1.x() < s2.x() : s1.y() < s2.y();
}

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
  QSize  size = binary.size();
  
  /* build a vector to hold images in Hough-space for radius 1..max_r, where
  max_r is specified or the maximum radius of a circle in this image */
  if(min_r == 0) min_r = 5;
  if(max_r == 0) max_r = MIN(source.width(), source.height()) / 2;
  
  QVector<Image> houghs(max_r - min_r);

  /* find all the edges */
  PointArray edge;
  for(unsigned int y = 0; y < size.height(); y++)
    for(unsigned int x = 0; x < size.width(); x++)
      if(binary.pixelIndex(x, y) == 1)
        edge.append(QPoint(x, y));

  /* construct vectorized boundary */
  {
      int low[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
      lower_bounds = (v8si *) low;

      int w = size.width();
      int h = size.height();
      int high[8] = { w, h, w, h, w, h, w, h };
      higher_bounds = (v8si *) high;
  }


  #pragma omp parallel for
  for(unsigned int i = min_r; i < max_r; i++)
  {
    /* instantiate circle template */
    const PointArray circle = circle_template(i);

    /* instantiate Hough-space for circles of radius i */
    Image &hough = houghs[i - min_r];
    hough.resize(size.width() * size.height());
    hough.fill(0);

    for (unsigned int k = 0; k < edge.size(); k++)
      accum_circle(hough, size, edge.at(k), circle);

    /* loop through all the Hough-space images, searching for bright spots, which
    indicate the center of a circle, then draw circles in image-space */
    unsigned int threshold = 4.9 * i;
    unsigned int total = size.width() * size.height();
    for (unsigned int k = 0; k < total; k++)
    {
      if (hough[k] > threshold)
      {
        unsigned int x = k / size.height();
        unsigned int y = k % size.height();
        draw_circle(detection, QPoint(x, y), circle, Qt::yellow);
      }
    }
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
void HoughCircleDetector::accum_circle(Image &image, const QSize &size, const QPoint &position, const PointArray &points)
{
  v8si *center, *offset; 
  v8si result;

  // set vectorized center
  PointArray centers(4);
  centers.fill(position);
  center = (v8si *) centers.constData();

  const int &w = size.width();
  const int &h = size.height();

  // iterate through circle edge points
  unsigned int total = points.size();
  for (int i = 0; i < total; i += 4) {
    offset = (v8si *) &points[i];
    result = *center + *offset;
#ifdef GCC_5_2
    v8si v1 = result > *lower_bounds;
    v8si v2 = result < *higher_bounds;
    v8si valid = v1 && v2;

    int *pos = (int *) &result;
    if (valid[0] && valid[1]) image[pos[0] * h + pos[1]]++;
    if (valid[2] && valid[3]) image[pos[2] * h + pos[3]]++;
    if (valid[4] && valid[5]) image[pos[4] * h + pos[5]]++;
    if (valid[6] && valid[7]) image[pos[6] * h + pos[7]]++;
#else
    v8si v1 = result - *lower_bounds;
    v8si v2 = result - *higher_bounds;
    int *pos = (int *) &result;

    if (v1[0] > 0 && v1[1] > 0 && v2[0] < 0 && v2[1] < 0) image[pos[0] * h + pos[1]]++;
    if (v1[2] > 0 && v1[3] > 0 && v2[2] < 0 && v2[3] < 0) image[pos[2] * h + pos[3]]++;
    if (v1[4] > 0 && v1[5] > 0 && v2[4] < 0 && v2[5] < 0) image[pos[4] * h + pos[5]]++;
    if (v1[6] > 0 && v1[7] > 0 && v2[6] < 0 && v2[7] < 0) image[pos[6] * h + pos[7]]++;
#endif
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
void HoughCircleDetector::draw_circle(QImage &image, const QPoint &position, const PointArray &points, const QColor &color)
{
  QRgb rgb = color.rgb();

  v8si *center, *offset; 
  v8si result;

  // set vectorized center
  PointArray centers(4);
  centers.fill(position);
  center = (v8si *) centers.constData();

  // iterate through circle edge points
  unsigned int total = points.size();
  for (int i = 0; i < total; i += 4) {
    offset = (v8si *) &points[i];
    result = *center + *offset;
#ifdef GCC_5_2
    v8si v1 = result > *lower_bounds;
    v8si v2 = result < *higher_bounds;
    v8si valid = v1 && v2;

    QPoint *pos = (QPoint *) &result;
    if (valid[0] && valid[1]) image.setPixel(pos[0], rgb);
    if (valid[2] && valid[3]) image.setPixel(pos[1], rgb);
    if (valid[4] && valid[5]) image.setPixel(pos[2], rgb);
    if (valid[6] && valid[7]) image.setPixel(pos[3], rgb);
#else
    v8si v1 = result - *lower_bounds;
    v8si v2 = result - *higher_bounds;

    QPoint *pos = (QPoint *) &result;
    if (v1[0] > 0 && v1[1] > 0 && v2[0] < 0 && v2[1] < 0) image.setPixel(pos[0], rgb);
    if (v1[2] > 0 && v1[3] > 0 && v2[2] < 0 && v2[3] < 0) image.setPixel(pos[1], rgb);
    if (v1[4] > 0 && v1[5] > 0 && v2[4] < 0 && v2[5] < 0) image.setPixel(pos[2], rgb);
    if (v1[6] > 0 && v1[7] > 0 && v2[6] < 0 && v2[7] < 0) image.setPixel(pos[3], rgb);
#endif
  }
}

const PointArray HoughCircleDetector::circle_template(unsigned int radius)
{
  int f = 1 - radius;
  int ddF_x = 1;
  int ddF_y = -2 * radius;
  int x = 0;
  int y = radius;

  PointArray points;
  points << QPoint(0,  radius)
         << QPoint(0, -radius)
         << QPoint( radius, 0)
         << QPoint(-radius, 0);

  while(x < y)
  {
    if(f >= 0)
    {
      y--;
      ddF_y += 2;
      f += ddF_y;
    }
    
    x++;
    ddF_x += 2;
    f += ddF_x;
    points << QPoint(+x, +y)
           << QPoint(+x, -y)
           << QPoint(-x, +y)
           << QPoint(-x, -y)
           << QPoint(+y, +x)
           << QPoint(+y, -x)
           << QPoint(-y, +x)
           << QPoint(-y, -x);
  }

  qSort(points.begin(), points.end(), QPointLessThan);
  return points;
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
  int w = source.width();
  int h = source.height();
  
  /*** Sobel edge detection ***/
  
  /* Set up Lx, Ly */
  QVector<QByteArray> Lx(3), Ly(3);
  
  Lx[0][0] = -1;  Lx[0][1] = +0;  Lx[0][2] = +1;
  Lx[1][0] = -2;  Lx[1][1] = +0;  Lx[1][2] = +2;
  Lx[2][0] = -1;  Lx[2][1] = +0;  Lx[2][2] = +1;
  
  Ly[0][0] = +1;  Ly[0][1] = +2;  Ly[0][2] = +1;
  Ly[1][0] = +0;  Ly[1][1] = +0;  Ly[1][2] = +0;
  Ly[2][0] = -1;  Ly[2][1] = -2;  Ly[2][2] = -1;
  
  #pragma omp parallel for
  for(unsigned int y = 0; y < h; y++)
  {
    #pragma omp parallel for
    for(unsigned int x = 0; x < w; x++)
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
          _x = _x < 0 ? -_x : _x >= w ? 2 * w - _x - 2 : _x;
          _y = _y < 0 ? -_y : _y >= h ? 2 * h - _y - 2 : _y;

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
  
  return binary;
}


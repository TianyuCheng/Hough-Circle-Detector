/****************************************************************************
**
** For Copyright & Licensing information, see COPYRIGHT in project root
**
****************************************************************************/

#include <QImage>
#include <QVector>

typedef QVector<unsigned int> IntArray;
typedef QVector<QPoint>       PointArray;
typedef QVector<IntArray>     Image;

class HoughCircleDetector
{
  public: /* class */
  
    HoughCircleDetector() {}
   ~HoughCircleDetector() {}
  
  public: /* methods */
  
    QImage detect(const QImage &source, unsigned int min_r, unsigned int max_r);
  
  private: /* methods */
  
    void accum_circle(Image &image, const QPoint &position, unsigned int radius);

    void draw_circle(QImage &image, const QPoint &position, unsigned int radius, const QColor &color);
    void draw_pixel(QImage &image, const int x, const int y, const QRgb &rgb);
    
    QImage edges(const QImage &source);
};


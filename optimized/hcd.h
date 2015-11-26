/****************************************************************************
**
** For Copyright & Licensing information, see COPYRIGHT in project root
**
****************************************************************************/

#include <QImage>
#include <QVector>

typedef QVector<unsigned int> IntArray;
typedef QVector<QPoint>       PointArray;
typedef IntArray              Image;

class HoughCircleDetector
{
  public: /* class */
  
    HoughCircleDetector() {}
   ~HoughCircleDetector() {}
  
  public: /* methods */
  
    QImage detect(const QImage &source, unsigned int min_r, unsigned int max_r);
  
  private: /* methods */
  
    void accum_circle(Image &image, const QSize &size, const QPoint &position, const PointArray &points);
    void draw_circle(QImage &image, const QPoint &position, const PointArray &points, const QColor &color);
    const PointArray circle_template(unsigned int radius);
    
    QImage edges(const QImage &source);
};


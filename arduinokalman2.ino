#include <math.h>

float quat[4] = {1.0f,0.0f,0.0f,0.0f};
float beta = 0.6f;
float zVariance = 500.0f;
float zAccelVariance = 1.0f;
float zAccelBiasVariance = 1.0f; 
float zInitial = 300; // comes from pressure sensor
float vInitial = 0.0f;
float aBiasInitial = 0.0f;

float z; 
float a;
float dt;  ////// elapsed time secs between measurements, can be calculated via function(while millis)
float *pZ;
float *pV;

float za;

float z_;  // position
float v_;  // velocity
float aBias_;  // acceleration

float zAccelBiasVariance_; // assumed fixed.
float zAccelVariance_;  // dynamic acceleration variance
float zVariance_; //  z measurement noise variance fixed

// constant values
float Pzz_ = 1.0f;
float Pzv_ = 0.0f;
float Pza_ = 0.0f;
float Pvz_ = 0.0f;
float Pvv_ = 1.0f;
float Pva_ = 0.0f;
float Paz_ = 0.0f;
float Pav_ = 0.0;
float Paa_ = 100000.0f;
float w = 0;
float ret = 0;
float scalefactor = 980.0f;

float fax;
float fay;
float faz;

float fmx;
float fmy;
float fmz;

float fgx;
float fgy;
float fgz;

float fa;
int bUseAccel;

float ax;
float ay;
float az;

float mx;
float my;
float mz;

float gx;
float gy;
float gz;

float q1 = quat[0];
float q2 = quat[1]; 
float q3 = quat[2]; 
float q4 = quat[3];

float norm;
float hx, hy, _2bx, _2bz;
float s1, s2, s3, s4;
float qDot1, qDot2, qDot3, qDot4;

float accel;
float _2q1mx;
float _2q1my;
float _2q1mz;
float _2q2mx;
float _4bx;
float _4bz;


float _2q1;
float _2q2;
float _2q3;
float _2q4;
float _2q1q3;
float _2q3q4;

float q1q1;
float q1q2;
float q1q3;
float q1q4;
float q2q2;
float q2q3;
float q2q4;
float q3q3;
float q3q4;
float q4q4;

float t00,t01,t02;
float t10,t11,t12;
float t20,t21,t22;
	
float dt2div2;
float dt3div2;
float dt4div4;

float innov;
float sInv;

float kz; 
float kv;
float ka;

void setup() {

    zAccelVariance_ = zAccelVariance;
    zAccelBiasVariance_ = zAccelBiasVariance;
    zVariance_ = zVariance;

    // values that given by user
    z_ = zInitial;
    v_ = vInitial;
    aBias_ = aBiasInitial;

}

void loop() {

    fa = sqrt(fax*fax + fay*fay + faz*faz); 
    bUseAccel = ((fa > 0.5f) && (fa < 1.5f)) ? 1 : 0;
    
/*  Modify according to your sensors
    fax = /////////////sensor read///////////////
    fay = /////////////sensor read//////////////
    faz = /////////////sensor read//////////////
    
    fmx = //////////////sensor read/////////////
    fmy = ///////////////sensor read////////////
    fmz = ///////////////sensor read////////////
    
    fgx = ///////////////sensor read////////////
    fgy = ///////////////sensor read////////////
    fgz = ///////////////sensor read///////////
*/

    ax = fax;
    ay = fay;
    az = faz;

    mx = fmx;
    my = fmy;
    mz = fmz;

    gx = fgx;  //*PI_DIV_180;
    gy = fgy;  //*PI_DIV_180;
    gz = fgz;  //*PI_DIV_180;

    // Compute rate of change of quaternion
    qDot1 = 0.5f * (-q2 * gx - q3 * gy - q4 * gz);
    qDot2 = 0.5f * (q1 * gx + q3 * gz - q4 * gy);
    qDot3 = 0.5f * (q1 * gy - q2 * gz + q4 * gx);
    qDot4 = 0.5f * (q1 * gz + q2 * gy - q3 * gx);

   if (bUseAccel) {

    _2q1 = 2.0f * q1;
    _2q2 = 2.0f * q2;
    _2q3 = 2.0f * q3;
    _2q4 = 2.0f * q4;
    _2q1q3 = 2.0f * q1 * q3;
    _2q3q4 = 2.0f * q3 * q4;

    q1q1 = q1 * q1;
    q1q2 = q1 * q2;
    q1q3 = q1 * q3;
    q1q4 = q1 * q4;
    q2q2 = q2 * q2;
    q2q3 = q2 * q3;
    q2q4 = q2 * q4;
    q3q3 = q3 * q3;
    q3q4 = q3 * q4;
    q4q4 = q4 * q4;


    // Normalise accelerometer measurement
    norm = sqrt(ax * ax + ay * ay + az * az);
    if (norm == 0.0f) return; // handle NaN
    norm = 1.0f/norm;
    ax *= norm;
    ay *= norm;
    az *= norm;

    // Normalise magnetometer measurement
    norm = sqrt(mx * mx + my * my + mz * mz);
    if (norm == 0.0f) return; // handle NaN
    norm = 1.0f/norm;
    mx *= norm;
    my *= norm;
    mz *= norm;

    // Reference direction of Earth's magnetic field
    _2q1mx = 2.0f * q1 * mx;
    _2q1my = 2.0f * q1 * my;
    _2q1mz = 2.0f * q1 * mz;
    _2q2mx = 2.0f * q2 * mx;
    hx = mx * q1q1 - _2q1my * q4 + _2q1mz * q3 + mx * q2q2 + _2q2 * my * q3 + _2q2 * mz * q4 - mx * q3q3 - mx * q4q4;
    hy = _2q1mx * q4 + my * q1q1 - _2q1mz * q2 + _2q2mx * q3 - my * q2q2 + my * q3q3 + _2q3 * mz * q4 - my * q4q4;
    _2bx = sqrt(hx * hx + hy * hy);
    _2bz = -_2q1mx * q3 + _2q1my * q2 + mz * q1q1 + _2q2mx * q4 - mz * q2q2 + _2q3 * my * q4 - mz * q3q3 + mz * q4q4;
    _4bx = 2.0f * _2bx;
    _4bz = 2.0f * _2bz;

    // Gradient descent algorithm corrective step
    s1 = -_2q3 * (2.0f * q2q4 - _2q1q3 - ax) + _2q2 * (2.0f * q1q2 + _2q3q4 - ay) - _2bz * q3 * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (-_2bx * q4 + _2bz * q2) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + _2bx * q3 * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
    s2 = _2q4 * (2.0f * q2q4 - _2q1q3 - ax) + _2q1 * (2.0f * q1q2 + _2q3q4 - ay) - 4.0f * q2 * (1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az) + _2bz * q4 * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (_2bx * q3 + _2bz * q1) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + (_2bx * q4 - _4bz * q2) * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
    s3 = -_2q1 * (2.0f * q2q4 - _2q1q3 - ax) + _2q4 * (2.0f * q1q2 + _2q3q4 - ay) - 4.0f * q3 * (1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az) + (-_4bx * q3 - _2bz * q1) * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (_2bx * q2 + _2bz * q4) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + (_2bx * q1 - _4bz * q3) * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
    s4 = _2q2 * (2.0f * q2q4 - _2q1q3 - ax) + _2q3 * (2.0f * q1q2 + _2q3q4 - ay) + (-_4bx * q4 + _2bz * q2) * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (-_2bx * q1 + _2bz * q3) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + _2bx * q2 * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
    norm = sqrt(s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4);    // normalise step magnitude
    norm = 1.0f/norm;
    s1 *= norm;
    s2 *= norm;
    s3 *= norm;
    s4 *= norm;

    // Compute rate of change of quaternion
    qDot1 -= beta * s1;
    qDot2 -= beta * s2;
    qDot3 -= beta * s3;
    qDot4 -= beta * s4;
   }

    // Integrate to yield quaternion
    q1 += qDot1 * dt;
    q2 += qDot2 * dt;
    q3 += qDot3 * dt;
    q4 += qDot4 * dt;
    norm = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);    // normalise quaternion
    norm = 1.0f/norm;
    quat[0] = q1 * norm;
    quat[1] = q2 * norm;
    quat[2] = q3 * norm;
    quat[3] = q4 * norm;

    za = 2.0f*(quat[1]*quat[3] - quat[0]*quat[2])*ax + 2.0f*(quat[0]*quat[1] + quat[2]*quat[3])*ay + (quat[0]*quat[0] - quat[1]*quat[1] - quat[2]*quat[2] + quat[3]*quat[3])*az - 1.0f;

    a = scalefactor*za;

    // Predict state
    accel = a - aBias_;
    v_ += accel * dt;
    z_ += v_ * dt;

    if (accel < 0) {
      w = -1 * accel;
    } else {
      w = accel;
    }

    zAccelVariance_ = w/50.0f;

    if (zAccelVariance_ < 1) {
      ret = 1;
    } else if (zAccelVariance_ > 50) {
      ret = 50;
    } else {
      ret = zAccelVariance_;
    }

    zAccelVariance_ = ret;

    // Predict State Covariance matrix
    dt2div2 = dt*dt/2.0f;
    dt3div2 = dt2div2*dt;
    dt4div4 = dt2div2*dt2div2;
	
    t00 = Pzz_ + dt*Pvz_ - dt2div2*Paz_;
    t01 = Pzv_ + dt*Pvv_ - dt2div2*Pav_;
    t02 = Pza_ + dt*Pva_ - dt2div2*Paa_;

    t10 = Pvz_ - dt*Paz_;
    t11 = Pvv_ - dt*Pav_;
    t12 = Pva_ - dt*Paa_;

    t20 = Paz_;
    t21 = Pav_;
    t22 = Paa_;

    Pzz_ = t00 + dt*t01 - dt2div2*t02;
    Pzv_ = t01 - dt*t02;
    Pza_ = t02;
	
    Pvz_ = t10 + dt*t11 - dt2div2*t12;
    Pvv_ = t11 - dt*t12;
    Pva_ = t12;
	
    Paz_ = t20 + dt*t21 - dt2div2*t22;
    Pav_ = t21 - dt*t22;
    Paa_ = t22;

    Pzz_ += dt4div4*zAccelVariance_;
    Pzv_ += dt3div2*zAccelVariance_;

    Pvz_ += dt3div2*zAccelVariance_;
    Pvv_ += dt*dt*zAccelVariance_;

    Paa_ += zAccelBiasVariance_;

    // Error
    innov = z - z_; 
    sInv = 1.0f / (Pzz_ + zVariance_);  

    // Kalman gains
    kz = Pzz_ * sInv;
    kv = Pvz_ * sInv;
    ka = Paz_ * sInv;

    // Update state 
    z_ += kz * innov;
    v_ += kv * innov;
    aBias_ += ka * innov;

    *pZ = z_;
    *pV = v_;

    // Update state covariance matrix
    Pzz_ -= kz * Pzz_;
    Pzv_ -= kz * Pzv_;
    Pza_ -= kz * Pza_;
	
    Pvz_ -= kv * Pzz_;
    Pvv_ -= kv * Pzv_;
    Pva_ -= kv * Pza_;
	
    Paz_ -= ka * Pzz_;
    Pav_ -= ka * Pzv_;
    Paa_ -= ka * Pza_;

}

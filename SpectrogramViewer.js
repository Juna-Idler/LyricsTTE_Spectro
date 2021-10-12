
class SpectrogramViewer
{
    constructor(audioBuffer,n = 1024 * 1,rate = 1.5,add = 1.0)
    {
        if (!Number.isInteger(Math.log2(n)))
            return
        this.N = n;
        if (!this.Window || this.Window.length != this.N)
        {
            this.Window = new Array(this.N);
            for (let i = 0; i < this.N;i++)
            {
                this.Window[i] = SpectrogramViewer.Hanning(i,this.N);
            }
        }
        if (!this.LogTable || this.LogTable.length != this.N)
        {
            this.LogTable = new Array(this.N / 2);
            this.CountTable = new Array(256);
            this.CountTable.fill(0);
            const max = Math.log10(this.N / 2);
            for (let i = 0; i < this.N / 2;i++)
            {
                this.LogTable[i] = Math.floor(Math.log10(i + 1) / max * 255);
                this.CountTable[this.LogTable[i]]++;
            }
            for (let i = 0; i < 256;i++)
            {
                if (this.CountTable[i] == 0)
                    this.CountTable[i] = 1;
            }
        }

        console.log("SpectrogramViewer start:duration=" + audioBuffer.duration + ",N=" +this.N);
        const start_time = performance.now();
        let mono;
        if (audioBuffer.numberOfChannels > 1)
        {
            mono = new Float32Array(audioBuffer.length);
            const l = audioBuffer.getChannelData(0);
            const r = audioBuffer.getChannelData(1);
            for (let i = 0; i < audioBuffer.length;i++)
            {
                mono[i] = (l[i] + r[i]) / 2;
            }
        }
        else if (audioBuffer.numberOfChannels > 0)
        {
            mono = audioBuffer.getChannelData(0);
        }

        this.width = Math.floor(audioBuffer.duration * 100);
        this.SpectrSet = new Array(this.width);
        
        const frequency = audioBuffer.sampleRate;

        let ar;
        let ai = new Float32Array(this.N);
        for (let i = 0,start = 0; start < mono.length - (frequency / 100) ;i++,start += frequency / 100)
        {
            ar = mono.slice(start,start + this.N);
            for (let j = 0;j < this.N;j++)
            {
                ar[j] *= this.Window[j];
            }
            ai.fill(0);
            SpectrogramViewer.FFT(ar,ai,this.N,false);
            const data = new Array(256);
            data.fill(0);
            for (let j = 0; j < this.N / 2;j++)
            {
                const v = Math.log10(Math.sqrt((ar[j] * ar[j] + ai[j] * ai[j]) / 2)) / rate + add;
                data[this.LogTable[j]] += v;
            }
            for (let j = 0; j < 256;j++)
            {
                data[j] /= this.CountTable[j];
            }
            this.SpectrSet[i] = data;
        }
        this.duration = performance.now() - start_time;
        console.log("SpectrogramViewer constructor end:" + this.duration);

    }

    get isValid() {return !!this.width;}

    DrawCanvas(canvas,start)
    {
        if (!this.isValid)
            return ;
        start = Math.floor(start);
        const ctx = canvas.getContext('2d');
       
        let width  = canvas.width;
        let height = 256;

//        ctx.fillStyle = "gray";
//        ctx.fillRect(0,0,width,height);
        ctx.clearRect(0,0,width,height);


        let x_start = 0;
        start = start | 0;
        if (start + width > this.width)
        {
            width = (this.width - start) | 0;
        }
        if (start < 0)
        {
            x_start = (-start) | 0;
            start = 0;
        }
        if (width <= 0)
            return;

        const imageData = ctx.createImageData(width,height);
        const data = imageData.data;
        data.fill(0);

        for (let x = x_start;x < width;x++)
        {
            for (let y = 0;y < 256;y++)
            {
                let v = this.SpectrSet[start + x - x_start][y];
                data[(((256 -y)*(imageData.width*4)) + (x*4)) + 0] = 255 * v;
                data[(((256 -y)*(imageData.width*4)) + (x*4)) + 1] = 255 * v;
                data[(((256 -y)*(imageData.width*4)) + (x*4)) + 2] = 255 * v;
                data[((y*(imageData.width*4)) + (x*4)) + 3] = 255;
            }
        }
        ctx.putImageData(imageData,0,0);
    }

    static FFT( an, bn, N, Inverse ){
        /////////////////////////////////////////
        //参考：Javaで学ぶシミュレーションの基礎
        /////////////////////////////////////////
        // 入力
        // N  ： 項数（2のべき乗）
        // an : 実数配列（順変換：実数空間データを項数Nで指定、逆変換：展開係数a(n)）
        // bn : 実数配列（順変換：虚数空間データを項数Nで指定、逆変換：展開係数b(n)）
        // Inverse : 逆変換の場合に true
        /////////////////////////////////////////
        // 出力
        // an : 実数配列（順変換：展開係数a(n)、逆変換：実数空間データ）
        // bn : 実数配列（順変換：展開係数b(n)、逆変換：虚数空間データ）
        /////////////////////////////////////////
        const ff = Inverse ? 1 : -1;
        const rot = new Array(N);
        for( let i = 0; i < rot.length; i++ ) rot[ i ] = 0;
        const nhalf = N/2;
        let num = N/2;
        const sc = 2 * Math.PI / N;
        while( num >= 1 ){
            for(let j = 0; j < N; j += 2 * num ){
                const phi = rot[j] / 2;
                const phi0 = phi + nhalf;
                const c = Math.cos( sc * phi );
                const s = Math.sin( sc * phi * ff );
                for( var k = j; k < j + num; k++ ){
                    const k1 = k + num;
                    const a0 = an[ k1 ] * c - bn[ k1 ] *s;
                    const b0 = an[ k1 ] * s + bn[ k1 ] *c;
                    an[ k1 ] = an[ k ] - a0;
                    bn[ k1 ] = bn[ k ] - b0;
                    an[ k ] = an[ k ] + a0;
                    bn[ k ] = bn[ k ] + b0;
                    rot[ k ] = phi;
                    rot[ k1 ] = phi0;
                }
            }
            num = num / 2;
        }
        for( let i = 0; i < N ; i++ ){
            const j = rot[ i ]; 
            if( j > i ){
                let tmp = an[ i ];
                an[ i ] = an[ j ];
                an[ j ] = tmp;
                tmp = bn[ i ];
                bn[ i ] = bn[ j ];
                bn[ j ] = tmp;  
            }
        }
        for( let i = 0; i < N ; i++ ){
            an[ i ] = an[ i ] / Math.sqrt(N);
            bn[ i ] = bn[ i ] / Math.sqrt(N);
        }
    }

    static Hanning(n,N)
    {
        return 0.5 - 0.5 * Math.cos( 2* Math.PI * n / (N-1) );
    }
}